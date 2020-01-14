import os
import sys
import argparse
import javabridge
import bioformats
import matplotlib.pyplot as plt
import glob
import itertools
from sklearn import svm
import skimage.filters
from scipy import ndimage as ndi
from skimage import restoration
from sklearn.cluster import KMeans
from skimage.feature import peak_local_max
from skimage import color
import pandas as pd
import numpy as np
import re
from skimage import exposure
from scipy.signal import argrelextrema
import joblib
from neighbor2d import line_profile_2d_v2
# from neighbor import line_profile_v2
from scipy.ndimage import binary_fill_holes
from scipy.ndimage import binary_opening
from matplotlib.colors import hsv_to_rgb
from skimage.future import graph
from ete3 import NCBITaxa
from skimage.filters import threshold_sauvola

javabridge.start_vm(class_path=bioformats.JARS)


def plot_taxon_color(taxon_lookup, input_folder):
    fig = plt.figure()
    fig.set_size_inches(12, 5)
    delta = 0.1
    for i in range(taxon_lookup.shape[0]):
        taxon_x, taxon_y = np.divmod(i, 10)
        plt.plot(taxon_x, taxon_y, 'o', color=tuple(
            hsv_to_rgb(taxon_lookup.loc[i, ['H', 'S', 'V']].values)))
        plt.text(taxon_x + delta, taxon_y - delta, taxon_lookup.loc[i, 'sci_name'], fontsize=8)
    plt.plot(taxon_x + 1, taxon_y, color=(1, 1, 1))
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('{}/taxon_color_lookup.png'.format(input_folder))
    plt.close()
    return


def convert_code_to_7b(code):
    bits = list(code)
    converted_code = ''.join([bits[i] for i in [0, 2, 3, 4, 7, 8, 9]])
    return(converted_code)


def convert_code_to_10b(code):
    bits = list(code)
    bits.insert(4, '0')
    bits.insert(4, '0')
    bits.insert(1, '0')
    converted_code = ''.join(bits)
    return(converted_code)


def load_ztslice(filename, z_index, t_index):
    image = bioformats.load_image(filename, z=z_index, t=t_index)
    return(image)


def get_t_range(filename):
    xml = bioformats.get_omexml_metadata(filename)
    ome = bioformats.OMEXML(xml)
    t_range = ome.image(0).Pixels.get_SizeT()
    return(t_range)


def get_z_range(filename):
    xml = bioformats.get_omexml_metadata(filename)
    ome = bioformats.OMEXML(xml)
    z_range = ome.image(0).Pixels.get_SizeZ()
    return(z_range)


def load_image_zstack_fixed_t(filename, t):
    z_range = get_z_range(filename)
    image = np.stack([load_ztslice(filename, k, t) for k in range(0, z_range)], axis=2)
    return(image)


def get_registered_image_from_tstack(filename):
    image_0 = load_image_zstack_fixed_t(filename, 0)
    image_registered = image_0.copy()
    image_0_sum = np.sum(image_0, axis=3)
    shift_vector_list = []
    nt = get_t_range(filename)
    for i in range(1, nt):
        image_i = load_image_zstack_fixed_t(filename, i)
        image_i_sum = np.sum(image_i, axis=3)
        shift_vector = skimage.feature.register_translation(image_0_sum, image_i_sum)[0]
        shift_vector = np.insert(shift_vector, 3, 0)
        shift_filter_mask = np.full(
            (image_0.shape[0], image_0.shape[1], image_0.shape[2]), False, dtype=bool)
        shift_x = int(shift_vector[0])
        shift_y = int(shift_vector[1])
        shift_z = int(shift_vector[2])
        original_x_min = int(np.maximum(0, shift_x))
        original_x_max = int(image_0.shape[0] + np.minimum(0, shift_x))
        original_y_min = int(np.maximum(0, shift_y))
        original_y_max = int(image_0.shape[1] + np.minimum(0, shift_y))
        original_z_min = int(np.maximum(0, shift_z))
        original_z_max = int(image_0.shape[2] + np.minimum(0, shift_z))
        registered_x_min = int(-np.minimum(0, shift_x))
        registered_x_max = int(image_0.shape[0] - np.maximum(0, shift_x))
        registered_y_min = int(-np.minimum(0, shift_y))
        registered_y_max = int(image_0.shape[1] - np.maximum(0, shift_y))
        registered_z_min = int(-np.minimum(0, shift_z))
        registered_z_max = int(image_0.shape[2] - np.maximum(0, shift_z))
        image_registered_hold = np.zeros(image_0.shape)
        image_registered_hold[original_x_min: original_x_max, original_y_min: original_y_max, original_z_min: original_z_max,
                              :] = image_i[registered_x_min: registered_x_max, registered_y_min: registered_y_max, registered_z_min: registered_z_max, :]
        shift_filter_mask[original_x_min: original_x_max,
                          original_y_min: original_y_max, original_z_min: original_z_max] = True
        image_registered += image_registered_hold
    return(image_registered)


def save_segmentation(segmentation, sample):
    seg_color = color.label2rgb(segmentation, bg_label=0, bg_color=(0, 0, 0))
    fig = plt.figure(frameon=False)
    fig.set_size_inches(5, 5)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    fig.add_axes(ax)
    ax.imshow(seg_color)
    segfilename = sample + '_seg.png'
    fig.savefig(segfilename, dpi=300)
    plt.close()
    np.save(sample + '_seg', segmentation)
    return


def save_identification(image_identification, sample):
    # seg_color = color.label2rgb(image_identification, bg_label = 0, bg_color = (0,0,0))
    fig = plt.figure(frameon=False)
    fig.set_size_inches(5, 5)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    fig.add_axes(ax)
    ax.imshow(image_identification)
    segfilename = sample + '_identification.png'
    fig.savefig(segfilename, dpi=300)
    plt.close()
    return


def save_sum_images(image_final, sample):
    fig = plt.figure(frameon=False)
    fig.set_size_inches(5, 5)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    fig.add_axes(ax)
    ax.imshow(image_final, cmap='jet')
    segfilename = sample + '_sum.png'
    fig.savefig(segfilename, dpi=300)
    plt.close()
    return


def save_enhanced_images(image_final, sample):
    fig = plt.figure(frameon=False)
    fig.set_size_inches(5, 5)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    fig.add_axes(ax)
    ax.imshow(image_final, cmap='jet')
    segfilename = sample + '_enhanced.png'
    fig.savefig(segfilename, dpi=300)
    plt.close()
    return


def generate_2d_segmentation(sample):
    excitations = ['488', '514', '561', '633']
    image_name = ['{}_{}.czi'.format(sample, x) for x in excitations]
    image_stack = [bioformats.load_image(filename) for filename in image_name]
    image_sum = [np.sum(image, axis=2) for image in image_stack]
    shift_vectors = [skimage.feature.register_translation(
        np.log(image_sum[0]+1e-4), np.log(image_sum[i]))[0] for i in range(1, 4)]
    shift_vectors.insert(0, np.asarray([0.0, 0.0]))
    image_registered = [np.zeros(image.shape) for image in image_stack]
    shift_filter_mask = [np.full((image.shape[0], image.shape[1]), False,
                                 dtype=bool) for image in image_stack]
    image_shape = image_stack[0].shape[0]
    for i in range(len(image_stack)):
        shift_row = int(shift_vectors[i][0])
        shift_col = int(shift_vectors[i][1])
        if np.abs(shift_row) > 15:
            shift_row = 0
        if np.abs(shift_col) > 15:
            shift_col = 0
        print(shift_row, shift_col)
        original_row_min = int(np.maximum(0, shift_row))
        original_row_max = int(image_shape + np.minimum(0, shift_row))
        original_col_min = int(np.maximum(0, shift_col))
        original_col_max = int(image_shape + np.minimum(0, shift_col))
        registered_row_min = int(-np.minimum(0, shift_row))
        registered_row_max = int(image_shape - np.maximum(0, shift_row))
        registered_col_min = int(-np.minimum(0, shift_col))
        registered_col_max = int(image_shape - np.maximum(0, shift_col))
        image_registered[i][original_row_min: original_row_max, original_col_min: original_col_max,
                            :] = image_stack[i][registered_row_min: registered_row_max, registered_col_min: registered_col_max, :]
        shift_filter_mask[i][original_row_min: original_row_max,
                             original_col_min: original_col_max] = True
    image_channel = np.dstack(image_registered)
    image_registered_sum = np.sum(image_channel, axis=2)
    image_registered_sum_norm = image_registered_sum/np.max(image_registered_sum)
    image_noise_variance = skimage.restoration.estimate_sigma(image_registered_sum_norm)
    image_registered_sum_nl = skimage.restoration.denoise_nl_means(image_registered_sum_norm, h=0.02)
    image_padded = skimage.util.pad(image_registered_sum_nl, 5, mode='edge')
    image_lp = line_profile_2d_v2(image_padded.astype(np.float64), 11, 9)
    image_lp = np.nan_to_num(image_lp)
    image_lp_min = np.min(image_lp, axis=3)
    image_lp_max = np.max(image_lp, axis=3)
    image_lp_max = image_lp_max - image_lp_min
    image_lp = image_lp - image_lp_min[:, :, :, None]
    image_lp_rel_norm = image_lp/image_lp_max[:, :, :, None]
    image_lp_rnc = image_lp_rel_norm[:, :, :, 5]
    image_lprns = np.average(image_lp_rnc, axis=2)
    image_lprn_lq = np.percentile(image_lp_rnc, 25, axis=2)
    image_lprn_uq = np.percentile(image_lp_rnc, 75, axis=2)
    image_lprn_qcv = np.zeros(image_lprn_uq.shape)
    image_lprn_qcv_pre = (image_lprn_uq - image_lprn_lq)/(image_lprn_uq + image_lprn_lq + 1e-8)
    image_lprn_qcv[image_lprn_uq > 0] = image_lprn_qcv_pre[image_lprn_uq > 0]
    image_final = image_lprns*(1-image_lprn_qcv)
    intensity_rough_seg = KMeans(n_clusters=2, random_state=0).fit_predict(
        image_final.reshape(np.prod(image_final.shape), 1)).reshape(image_final.shape)
    image0 = image_final*(intensity_rough_seg == 0)
    image1 = image_final*(intensity_rough_seg == 1)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    if (i0 < i1):
        intensity_rough_seg_mask = (intensity_rough_seg == 1)
        intensity_rough_seg_bkg = (intensity_rough_seg == 0)
    else:
        intensity_rough_seg_mask = (intensity_rough_seg == 0)
        intensity_rough_seg_bkg = (intensity_rough_seg == 1)
    image_lprns_rsfbo = skimage.morphology.binary_opening(intensity_rough_seg_mask)
    image_lprns_rsfbosm = skimage.morphology.remove_small_objects(image_lprns_rsfbo, 10)
    image_lprns_rsfbosm_bfh = binary_fill_holes(image_lprns_rsfbosm)
    intensity_rough_seg_mask_bfh = binary_fill_holes(intensity_rough_seg_mask)
    image_watershed_seeds = skimage.measure.label(
        image_lprns_rsfbosm_bfh*intensity_rough_seg_mask_bfh)
    image_registered_sum_nl_log = np.log10(image_registered_sum_nl)
    image_bkg_filter = KMeans(n_clusters=2, random_state=0).fit_predict(image_registered_sum_nl_log.reshape(
        np.prod(image_registered_sum_nl_log.shape), 1)).reshape(image_registered_sum_nl_log.shape)
    image0 = image_registered_sum_nl*(image_bkg_filter == 0)
    image1 = image_registered_sum_nl*(image_bkg_filter == 1)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    if (i0 < i1):
        image_bkg_filter_mask = (image_bkg_filter == 1)
    else:
        image_bkg_filter_mask = (image_bkg_filter == 0)
    image_final_bkg_filtered = image_registered_sum_nl_log*image_bkg_filter_mask
    image_sum_bkg_filtered = image_registered_sum*image_bkg_filter_mask
    image_watershed_seeds_bkg_filtered = image_watershed_seeds*image_bkg_filter_mask
    image_watershed_sauvola = threshold_sauvola(image_registered_sum_nl)
    image_watershed_mask_bkg_filtered = (
        image_registered_sum_nl > image_watershed_sauvola)*image_bkg_filter_mask
    image_seg = skimage.morphology.watershed(-image_final_bkg_filtered,
                                             image_watershed_seeds_bkg_filtered, mask=image_watershed_mask_bkg_filtered)
    image_seg = skimage.segmentation.relabel_sequential(image_seg)[0]
    adjacency_seg = skimage.morphology.watershed(-image_sum_bkg_filtered,
                                                 image_watershed_seeds_bkg_filtered, mask=image_bkg_filter_mask)
    adjacency_seg = skimage.segmentation.relabel_sequential(adjacency_seg)[0]
    save_segmentation(image_seg, sample)
    return(image_registered_sum, image_channel, image_final_bkg_filtered, image_seg, adjacency_seg)


def generate_3d_segmentation(sample):
    excitations = ['488', '514', '561', '633']
    image_name = ['{}_{}.czi'.format(sample, x) for x in excitations]
    image_stack = [get_registered_image_from_tstack(filename) for filename in image_name]
    image_sum = [np.sum(image, axis=3) for image in image_stack]
    shift_vectors = [skimage.feature.register_translation(
        np.log(image_sum[0]), np.log(image_sum[i]))[0] for i in range(1, 4)]
    shift_vectors.insert(0, np.asarray([0.0, 0.0, 0.0]))
    image_registered = [np.zeros(image.shape) for image in image_stack]
    shift_filter_mask = [np.full((image.shape[0], image.shape[1], image.shape[2]),
                                 False, dtype=bool) for image in image_stack]
    image_shape = image_stack[0].shape
    for i in range(len(image_stack)):
        shift_x = int(shift_vectors[i][0])
        shift_y = int(shift_vectors[i][1])
        shift_z = int(shift_vectors[i][2])
        print(shift_x, shift_y, shift_z)
        original_x_min = int(np.maximum(0, shift_x))
        original_x_max = int(image_shape[0] + np.minimum(0, shift_x))
        original_y_min = int(np.maximum(0, shift_y))
        original_y_max = int(image_shape[1] + np.minimum(0, shift_y))
        original_z_min = int(np.maximum(0, shift_z))
        original_z_max = int(image_shape[2] + np.minimum(0, shift_z))
        registered_x_min = int(-np.minimum(0, shift_x))
        registered_x_max = int(image_shape[0] - np.maximum(0, shift_x))
        registered_y_min = int(-np.minimum(0, shift_y))
        registered_y_max = int(image_shape[1] - np.maximum(0, shift_y))
        registered_z_min = int(-np.minimum(0, shift_z))
        registered_z_max = int(image_shape[2] - np.maximum(0, shift_z))
        image_registered[i][original_x_min: original_x_max, original_y_min: original_y_max, original_z_min: original_z_max,
                            :] = image_stack[i][registered_x_min: registered_x_max, registered_y_min: registered_y_max, registered_z_min: registered_z_max, :]
        shift_filter_mask[i][original_x_min: original_x_max,
                             original_y_min: original_y_max, original_z_min: original_z_max] = True
    image_channel = np.concatenate(image_registered, axis=3)
    image_registered_sum = np.sum(image_channel, axis=3)
    image_registered_sum_norm = image_registered_sum/np.max(image_registered_sum)
    image_noise_variance = skimage.restoration.estimate_sigma(image_registered_sum_norm)
    image_registered_sum_nl = skimage.restoration.denoise_nl_means(image_registered_sum_norm, h=0.02)
    image_padded = skimage.util.pad(image_registered_sum_nl, 5, mode='edge')
    image_lp = line_profile_v2(image_padded.astype(np.float64), 11, 9, 9)
    image_lp = np.nan_to_num(image_lp)
    image_lp_min = np.min(image_lp, axis=4)
    image_lp_max = np.max(image_lp, axis=4)
    image_lp_max = image_lp_max - image_lp_min
    image_lp = image_lp - image_lp_min[:, :, :, :, None]
    image_lp_rel_norm = image_lp/image_lp_max[:, :, :, :, None]
    image_lp_rnc = image_lp_rel_norm[:, :, :, :, 5]
    image_lprns = np.average(image_lp_rnc, axis=3)
    image_lprn_lq = np.percentile(image_lp_rnc, 25, axis=3)
    image_lprn_uq = np.percentile(image_lp_rnc, 75, axis=3)
    image_lprn_qcv = (image_lprn_uq - image_lprn_lq)/(image_lprn_uq + image_lprn_lq)
    image_lprn_qcv = np.nan_to_num(image_lprn_qcv)
    image_final = image_lprns*(1-image_lprn_qcv)
    intensity_rough_seg = np.zeros(image_final.shape).astype(int)
    intensity_rough_seg[image_final > 0] = KMeans(n_clusters=2, random_state=0).fit_predict(
        image_final[image_final > 0].reshape(-1, 1))
    image0 = image_final[image_final > 0]*(intensity_rough_seg[image_final > 0] == 0)
    image1 = image_final[image_final > 0]*(intensity_rough_seg[image_final > 0] == 1)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    intensity_rough_seg_mask = np.zeros(image_final.shape).astype(int)
    if i0 < i1:
        intensity_rough_seg_mask[image_final > 0] = intensity_rough_seg[image_final > 0]
    else:
        intensity_rough_seg_mask[image_final > 0] = intensity_rough_seg[image_final > 0] == 0
    image_final_seg = np.zeros(image_final.shape).astype(int)
    image_final_seg[image_final > 0] = KMeans(n_clusters=3, random_state=0).fit_predict(
        image_final[image_final > 0].reshape(-1, 1))
    image0 = image_final[image_final > 0]*(image_final_seg[image_final > 0] == 0)
    image1 = image_final[image_final > 0]*(image_final_seg[image_final > 0] == 1)
    image2 = image_final[image_final > 0]*(image_final_seg[image_final > 0] == 2)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    i2 = np.average(image2[image2 > 0])
    image_lprns_rsf = np.zeros(image_final.shape).astype(int)
    image_lprns_rsf[image_final > 0] = (
        image_final_seg[image_final > 0] == np.argmax([i0, i1, i2]))*1
    image_lprns_rsfbo = binary_opening(image_lprns_rsf)
    image_lprns_rsfbosm = skimage.morphology.remove_small_objects(image_lprns_rsfbo, 10)
    image_lprns_rsfbosm_bfh = binary_fill_holes(image_lprns_rsfbosm)
    intensity_rough_seg_mask_bfh = binary_fill_holes(intensity_rough_seg_mask)
    image_watershed_mask = image_lprns_rsfbosm_bfh*intensity_rough_seg_mask_bfh
    cell_sm_label = skimage.morphology.label(image_watershed_mask)
    dist_lab = skimage.morphology.label(cell_sm_label)
    markers = skimage.measure.regionprops(dist_lab)
    dist_be = np.zeros(dist_lab.shape)
    while(len(markers) > 0):
        for j in range(0, len(markers)):
            a = markers[j].area
            if (a < 600):
                dist_be[dist_lab == j+1] = 1
                dist_lab[dist_lab == j+1] = 0
        dist_bin_temp = skimage.morphology.binary_erosion(dist_lab)
        dist_bin_temp_sm = skimage.morphology.remove_small_objects(dist_bin_temp, 10)
        dist_lab = skimage.morphology.label(dist_bin_temp_sm)
        markers = skimage.measure.regionprops(dist_lab)
    dist_final = skimage.morphology.label(
        skimage.morphology.remove_small_objects(skimage.morphology.label(dist_be), 10))
    image_watershed_seeds = skimage.morphology.label(dist_final)
    image_registered_sum_nl_log = np.log10(image_registered_sum_nl)
    image_bkg_filter = KMeans(n_clusters=2, random_state=0).fit_predict(image_registered_sum_nl_log.reshape(
        np.prod(image_registered_sum_nl_log.shape), 1)).reshape(image_registered_sum_nl_log.shape)
    image0 = image_registered_sum_nl*(image_bkg_filter == 0)
    image1 = image_registered_sum_nl*(image_bkg_filter == 1)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    if (i0 < i1):
        image_bkg_filter_mask = (image_bkg_filter == 1)
    else:
        image_bkg_filter_mask = (image_bkg_filter == 0)
    image_final_bkg_filtered = image_final*image_bkg_filter_mask
    image_watershed_seeds_bkg_filtered = image_watershed_seeds*image_bkg_filter_mask
    image_watershed_mask_bkg_filtered = image_watershed_mask*image_bkg_filter_mask
    image_seg = skimage.morphology.watershed(-image_final_bkg_filtered,
                                             image_watershed_seeds_bkg_filtered, mask=image_watershed_mask_bkg_filtered)
    image_seg = skimage.segmentation.relabel_sequential(image_seg)[0]
    # save_segmentation(image_seg, sample)
    return(image_registered_sum, image_channel, image_final_bkg_filtered, image_seg)


def generate_3d_segmentation_skeleton(sample):
    excitations = ['488', '514', '561', '633']
    image_name = ['{}_{}.czi'.format(sample, x) for x in excitations]
    image_stack = [get_registered_image_from_tstack(filename) for filename in image_name]
    image_sum = [np.sum(image, axis=3) for image in image_stack]
    shift_vectors = [skimage.feature.register_translation(
        np.log(image_sum[0]), np.log(image_sum[i]))[0] for i in range(1, 4)]
    shift_vectors.insert(0, np.asarray([0.0, 0.0, 0.0]))
    image_registered = [np.zeros(image.shape) for image in image_stack]
    shift_filter_mask = [np.full((image.shape[0], image.shape[1], image.shape[2]),
                                 False, dtype=bool) for image in image_stack]
    image_shape = image_stack[0].shape
    for i in range(len(image_stack)):
        shift_x = int(shift_vectors[i][0])
        shift_y = int(shift_vectors[i][1])
        shift_z = int(shift_vectors[i][2])
        print(shift_x, shift_y, shift_z)
        original_x_min = int(np.maximum(0, shift_x))
        original_x_max = int(image_shape[0] + np.minimum(0, shift_x))
        original_y_min = int(np.maximum(0, shift_y))
        original_y_max = int(image_shape[1] + np.minimum(0, shift_y))
        original_z_min = int(np.maximum(0, shift_z))
        original_z_max = int(image_shape[2] + np.minimum(0, shift_z))
        registered_x_min = int(-np.minimum(0, shift_x))
        registered_x_max = int(image_shape[0] - np.maximum(0, shift_x))
        registered_y_min = int(-np.minimum(0, shift_y))
        registered_y_max = int(image_shape[1] - np.maximum(0, shift_y))
        registered_z_min = int(-np.minimum(0, shift_z))
        registered_z_max = int(image_shape[2] - np.maximum(0, shift_z))
        image_registered[i][original_x_min: original_x_max, original_y_min: original_y_max, original_z_min: original_z_max,
                            :] = image_stack[i][registered_x_min: registered_x_max, registered_y_min: registered_y_max, registered_z_min: registered_z_max, :]
        shift_filter_mask[i][original_x_min: original_x_max,
                             original_y_min: original_y_max, original_z_min: original_z_max] = True
    image_channel = np.concatenate(image_registered, axis=3)
    image_registered_sum = np.sum(image_channel, axis=3)
    image_registered_sum_norm = image_registered_sum/np.max(image_registered_sum)
    image_noise_variance = skimage.restoration.estimate_sigma(image_registered_sum_norm)
    image_registered_sum_nl = skimage.restoration.denoise_nl_means(image_registered_sum_norm, h=0.02)
    image_padded = skimage.util.pad(image_registered_sum_nl, 5, mode='edge')
    image_lp = line_profile_v2(image_padded.astype(np.float64), 11, 9, 9)
    image_lp = np.nan_to_num(image_lp)
    image_lp_min = np.min(image_lp, axis=4)
    image_lp_max = np.max(image_lp, axis=4)
    image_lp_max = image_lp_max - image_lp_min
    image_lp = image_lp - image_lp_min[:, :, :, :, None]
    image_lp_rel_norm = image_lp/image_lp_max[:, :, :, :, None]
    image_lp_rnc = image_lp_rel_norm[:, :, :, :, 5]
    image_lprns = np.average(image_lp_rnc, axis=3)
    image_lprn_lq = np.percentile(image_lp_rnc, 25, axis=3)
    image_lprn_uq = np.percentile(image_lp_rnc, 75, axis=3)
    image_lprn_qcv = (image_lprn_uq - image_lprn_lq)/(image_lprn_uq + image_lprn_lq)
    image_lprn_qcv = np.nan_to_num(image_lprn_qcv)
    image_final = image_lprns*(1-image_lprn_qcv)
    intensity_rough_seg = np.zeros(image_final.shape).astype(int)
    intensity_rough_seg[image_final > 0] = KMeans(n_clusters=2, random_state=0).fit_predict(
        image_final[image_final > 0].reshape(-1, 1))
    image0 = image_final[image_final > 0]*(intensity_rough_seg[image_final > 0] == 0)
    image1 = image_final[image_final > 0]*(intensity_rough_seg[image_final > 0] == 1)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    intensity_rough_seg_mask = np.zeros(image_final.shape).astype(int)
    if i0 < i1:
        intensity_rough_seg_mask[image_final > 0] = intensity_rough_seg[image_final > 0]
    else:
        intensity_rough_seg_mask[image_final > 0] = intensity_rough_seg[image_final > 0] == 0
    image_final_seg = np.zeros(image_final.shape).astype(int)
    image_final_seg[image_final > 0] = KMeans(n_clusters=3, random_state=0).fit_predict(
        image_final[image_final > 0].reshape(-1, 1))
    image0 = image_final[image_final > 0]*(image_final_seg[image_final > 0] == 0)
    image1 = image_final[image_final > 0]*(image_final_seg[image_final > 0] == 1)
    image2 = image_final[image_final > 0]*(image_final_seg[image_final > 0] == 2)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    i2 = np.average(image2[image2 > 0])
    image_lprns_rsf = np.zeros(image_final.shape).astype(int)
    image_lprns_rsf[image_final > 0] = (
        image_final_seg[image_final > 0] == np.argmax([i0, i1, i2]))*1
    image_lprns_rsfbo = binary_opening(image_lprns_rsf)
    image_lprns_rsfbosm = skimage.morphology.remove_small_objects(image_lprns_rsfbo, 10)
    image_lprns_rsfbosm_bfh = binary_fill_holes(image_lprns_rsfbosm)
    intensity_rough_seg_mask_bfh = binary_fill_holes(intensity_rough_seg_mask)
    image_watershed_mask = image_lprns_rsfbosm_bfh*intensity_rough_seg_mask_bfh

    cell_sm_label = skimage.morphology.label(image_watershed_mask)
    dist_lab = skimage.morphology.label(cell_sm_label)
    markers = skimage.measure.regionprops(dist_lab)
    dist_be = np.zeros(dist_lab.shape)

    for j in range(0, len(markers)):
        print(j)
        marker_skeleton = skimage.morphology.skeletonize_3d(cell_sm_label == markers[j].label)
        marker_skeleton_padded = skimage.util.pad(marker_skeleton, 1, mode='constant')
        bx_min, by_min, bz_min, bx_max, by_max, bz_max = markers[j].bbox
        cell_sm_label_local = cell_sm_label[bx_min:bx_max, by_min:by_max, bz_min:bz_max]
        connectivity = np.zeros(cell_sm_label_local.shape)
        for i in range(cell_sm_label_local.shape[0]):
            for j in range(cell_sm_label_local.shape[1]):
                for k in range(cell_sm_label_local.shape[2]):
                    connectivity[i, j, k] = np.sum(
                        marker_skeleton_padded[bx_min+i-1:bx_min+i+2, by_min+j-1:by_min+j+2, bz_min+k-1:bz_min+k+2])
                    dist_be[bx_min:bx_max, by_min:by_max, bz_min:bz_max] = (connectivity == 3)*1

    image_watershed_seeds = skimage.morphology.label(
        skimage.morphology.remove_small_objects(dist_be, 10))
    image_registered_sum_nl_log = np.log10(image_registered_sum_nl)
    image_bkg_filter = KMeans(n_clusters=2, random_state=0).fit_predict(image_registered_sum_nl_log.reshape(
        np.prod(image_registered_sum_nl_log.shape), 1)).reshape(image_registered_sum_nl_log.shape)
    image0 = image_registered_sum_nl*(image_bkg_filter == 0)
    image1 = image_registered_sum_nl*(image_bkg_filter == 1)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    if (i0 < i1):
        image_bkg_filter_mask = (image_bkg_filter == 1)
    else:
        image_bkg_filter_mask = (image_bkg_filter == 0)
    image_final_bkg_filtered = image_final*image_bkg_filter_mask
    image_watershed_seeds_bkg_filtered = image_watershed_seeds*image_bkg_filter_mask
    image_watershed_mask_bkg_filtered = image_watershed_mask*image_bkg_filter_mask
    image_seg = skimage.morphology.watershed(-image_final_bkg_filtered,
                                             image_watershed_seeds_bkg_filtered, mask=image_watershed_mask_bkg_filtered)
    image_seg = skimage.segmentation.relabel_sequential(image_seg)[0]
    # save_segmentation(image_seg, sample)
    return(image_registered_sum, image_channel, image_final_bkg_filtered, image_seg)


def measure_biofilm_images_no_reference(sample, dimension):
    if dimension == 2:
        image_registered_sum, image_registered, image_final_bkg_filtered, segmentation, adjacency_seg = generate_2d_segmentation(
            sample)
        save_sum_images(image_registered_sum, sample)
        save_enhanced_images(image_final_bkg_filtered, sample)
    else:
        image_registered_sum, image_registered, image_final_bkg_filtered, segmentation = generate_3d_segmentation(
            sample)
    np.save('{}_seg.npy'.format(sample), segmentation)
    np.save('{}_registered.npy'.format(sample), image_registered)
    cells = skimage.measure.regionprops(segmentation)
    if dimension == 2:
        avgint = np.empty((len(cells), image_registered.shape[2]))
        for k in range(0, image_registered.shape[2]):
            cells = skimage.measure.regionprops(
                segmentation, intensity_image=image_registered[:, :, k])
            avgint[:, k] = [x.mean_intensity for x in cells]
    else:
        avgint = np.empty((len(cells), image_registered.shape[3]))
        for k in range(0, image_registered.shape[3]):
            cells = skimage.measure.regionprops(
                segmentation, intensity_image=image_registered[:, :, :, k])
            avgint[:, k] = [x.mean_intensity for x in cells]
    avgint_norm = avgint/np.max(avgint, axis=1)[:, None]
    pd.DataFrame(avgint_norm).to_csv('{}_avgint_norm.csv'.format(sample), index=None)
    pd.DataFrame(avgint).to_csv('{}_avgint.csv'.format(sample), index=None)
    return


def measure_biofilm_images(sample, dimension, umap_transform, clf_umap, clf, taxon_lookup):
    if dimension == 2:
        image_registered_sum, image_registered, image_final_bkg_filtered, segmentation, adjacency_seg = generate_2d_segmentation(
            sample)
        save_sum_images(image_registered_sum, sample)
        save_enhanced_images(image_final_bkg_filtered, sample)
    else:
        image_registered_sum, image_registered, image_final_bkg_filtered, segmentation = generate_3d_segmentation(
            sample)
    np.save('{}_seg.npy'.format(sample), segmentation)
    np.save('{}_registered.npy'.format(sample), image_registered)
    cells = skimage.measure.regionprops(segmentation)
    if dimension == 2:
        avgint = np.empty((len(cells), image_registered.shape[2]))
        for k in range(0, image_registered.shape[2]):
            cells = skimage.measure.regionprops(
                segmentation, intensity_image=image_registered[:, :, k])
            avgint[:, k] = [x.mean_intensity for x in cells]
    else:
        avgint = np.empty((len(cells), image_registered.shape[3]))
        for k in range(0, image_registered.shape[3]):
            cells = skimage.measure.regionprops(
                segmentation, intensity_image=image_registered[:, :, :, k])
            avgint[:, k] = [x.mean_intensity for x in cells]
    avgint_norm = avgint/np.max(avgint, axis=1)[:, None]
    avgint_norm = np.concatenate((avgint_norm, np.zeros((avgint_norm.shape[0], 4))), axis=1)
    avgint_norm[:, 63] = clf[0].predict(avgint_norm[:, 0:23])
    avgint_norm[:, 64] = clf[1].predict(avgint_norm[:, 23:43])
    avgint_norm[:, 65] = clf[2].predict(avgint_norm[:, 43:57])
    avgint_norm[:, 66] = clf[3].predict(avgint_norm[:, 57:63])
    avgint_umap_transformed = umap_transform.transform(avgint_norm)
    cell_ids_norm = clf_umap.predict(avgint_umap_transformed)
    cell_info = pd.DataFrame(np.concatenate((avgint_norm, cell_ids_norm[:, None]), axis=1))
    cell_info[68] = sample
    cell_info[69] = np.asarray([x.label for x in cells])
    if dimension == 2:
        cell_info[70] = np.asarray([x.centroid[0] for x in cells])
        cell_info[71] = np.asarray([x.centroid[1] for x in cells])
        cell_info[72] = np.asarray([x.major_axis_length for x in cells])
        cell_info[73] = np.asarray([x.minor_axis_length for x in cells])
        cell_info[74] = np.asarray([x.eccentricity for x in cells])
        cell_info[75] = np.asarray([x.orientation for x in cells])
        cell_info[76] = np.asarray([x.area for x in cells])
        ids = list(set(cell_ids_norm))
        # ids_converted = [convert_code_to_7b(x) for x in ids]
        image_identification = np.zeros((segmentation.shape[0], segmentation.shape[1], 3))
        image_identification_barcode = np.zeros(segmentation.shape)
        for q in range(0, len(ids)):
            cell_population = np.where(cell_ids_norm == ids[q])[0]
            for r in range(0, len(cell_population)):
                image_identification_barcode[segmentation == cell_population[r]+1] = int(ids[q], 2)
                if ids[q] in taxon_lookup.code.values:
                    image_identification[segmentation == cell_population[r]+1, :] = hsv_to_rgb(
                        taxon_lookup.loc[taxon_lookup.code == ids[q], ['H', 'S', 'V']].values)
                else:
                    image_identification[segmentation ==
                                         cell_population[r]+1, :] = np.array([0, 0, 0])
        save_identification(image_identification, sample)
    else:
        cell_info[70] = np.asarray([x.centroid[0] for x in cells])
        cell_info[71] = np.asarray([x.centroid[1] for x in cells])
        cell_info[72] = np.asarray([x.centroid[2] for x in cells])
    cellinfofilename = sample + '_cell_information.csv'
    cell_info.to_csv(cellinfofilename, index=None, header=None)
    # edge_map = skimage.filters.sobel(segmentation > 0)
    # rag = skimage.future.graph.rag_boundary(adjacency_seg, edge_map)
    # adjacency_matrix = pd.DataFrame(np.zeros((taxon_lookup.shape[0], taxon_lookup.shape[0])), index = taxon_lookup.code.values, columns = taxon_lookup.code.values)
    # for i in range(cell_info.shape[0]):
    #     edges = list(rag.edges(i+1))
    #     for e in edges:
    #         node_1 = e[0]
    #         node_2 = e[1]
    #         if (node_1 != 0) and (node_2 !=0):
    #             barcode_1 = convert_code_to_7b(cell_info.iloc[node_1-1,67])
    #             barcode_2 = convert_code_to_7b(cell_info.iloc[node_2-1, 67])
    #             adjacency_matrix.loc[barcode_1, barcode_2] += 1
    # adjacencyfilename = sample + '_adjacency_matrix.csv'
    # adjacency_matrix.to_csv(adjacencyfilename)
    return


def main():
    parser = argparse.ArgumentParser('Design FISH probes for a complex microbial community')
    parser.add_argument('input_folder', type=str,
                        help='Input folder containing images of biological samples')
    parser.add_argument('-p', '--probe_design_filename', dest='probe_design_filename',
                        type=str, default='', help='Input folder containing images of biological samples')
    parser.add_argument('-r', '--ref_clf', dest='ref_clf', type=str, default='',
                        help='Input folder containing images of biological samples')
    parser.add_argument('-d', '--dimension', dest='dimension',
                        type=int, default='', help='Reference folder')
    parser.add_argument('-s', '--subfolder', dest='subfolder',
                        type=str, default='F', help='Sub folder')
    parser.add_argument('-e', '--epithelial', dest='ep', type=str, default='F', help='Sub folder')
    args = parser.parse_args()
    if args.probe_design_filename == '':
        filenames = glob.glob('{}/*.czi'.format(args.input_folder))
        samples = list(set([re.sub('_[0-9][0-9][0-9].czi', '', file) for file in filenames]))
        i = 1
        for s in samples:
            measure_biofilm_images_no_reference(s, args.dimension)
            print("Finished str(i) of str(len(samples))")
            i = i + 1
    else:
        probes = pd.read_csv(args.probe_design_filename, dtype={'code': str})
        ncbi = NCBITaxa()
        taxon_lookup = probes.loc[:, ['target_taxon', 'code']].drop_duplicates()
        taxon_lookup['H'] = np.arange(0, 1, 1/taxon_lookup.shape[0])
        taxon_lookup['S'] = 1
        taxon_lookup['V'] = 1
        taxon_sciname = pd.DataFrame.from_dict(ncbi.get_taxid_translator(
            taxon_lookup.target_taxon.values), orient='index').reset_index()
        taxon_sciname.columns = ['target_taxon', 'sci_name']
        taxon_lookup = taxon_lookup.merge(taxon_sciname, on='target_taxon')
        taxon_lookup.to_csv('{}/taxon_color_lookup.csv'.format(args.input_folder))
        if args.ep == 'T':
            taxon_lookup.loc[taxon_lookup.shape[0]] = ['0', '0000000', 0, 0, 0.5, 'Epithelial']
        umap_transform = joblib.load(args.ref_clf)
        clf_umap = joblib.load(re.sub('transform_biofilm_7b.pkl',
                                      'transformed_biofilm_7b_svc.pkl', args.ref_clf))
        clf = joblib.load(re.sub('transform_biofilm_7b.pkl',
                                 'transformed_biofilm_7b_check_svc.pkl', args.ref_clf))
        if args.subfolder == 'T':
            sf = glob.glob('{}/*'.format(args.input_folder))
            for subf in sf:
                filenames = glob.glob('{}/*.czi'.format(subf))
                samples = list(set([re.sub('_[0-9][0-9][0-9].czi', '', file) for file in filenames]))
                for s in samples:
                    measure_biofilm_images(s, args.dimension, umap_transform,
                                           clf_umap, clf, taxon_lookup)
        else:
            filenames = glob.glob('{}/*.czi'.format(args.input_folder))
            samples = list(set([re.sub('_[0-9][0-9][0-9].czi', '', file) for file in filenames]))
            for s in samples:
                measure_biofilm_images(s, args.dimension, umap_transform,
                                       clf_umap, clf, taxon_lookup)
    return


if __name__ == '__main__':
    main()

javabridge.kill_vm()

# skeleton_connectivity = np.zeros((500,500,50))
# skeleton = skeleton/np.max(skeleton)
#
# for i in range(500):
#     for j in range(500):
#         for k in range(50):
#             if skeleton[i,j,k] > 0:
#                 skeleton_connectivity[i,j,k] = np.sum(skeleton[i-1:i+2,j-1:j+2,k-1:k+2])

import javabridge
import bioformats
import numpy as np
from skimage.feature import register_translation

javabridge.start_vm(class_path=bioformats.JARS)

sample = 'data/images/2019_11_06_probe_DNA_saber_copynumber_high_probe_all_fov_1_chan'
excitations = ['488', '514', '561', '633']
image_name = ['{}_{}.czi'.format(sample, x) for x in excitations]

image_stack = [bioformats.load_image(filename) for filename in image_name]

image_sum = [np.sum(image, axis=2) for image in image_stack]
shift_vectors = [register_translation(np.log(image_sum[0]+1e-4), np.log(image_sum[i]))[0] for i in range(1, 4)]
shift_vectors.insert(0, np.asarray([0.0, 0.0, 0.0]))
image_registered = [np.zeros(image.shape) for image in image_stack]
shift_filter_mask = [np.full((image.shape[0], image.shape[1], image.shape[2]), False, dtype=bool) for image in image_stack]
image_shape = image_stack[0].shape
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
    image_registered[i][original_row_min: original_row_max, original_col_min: original_col_max, :] = image_stack[i][registered_row_min: registered_row_max, registered_col_min: registered_col_max, :]
    shift_filter_mask[i][original_row_min: original_row_max, original_col_min: original_col_max] = True

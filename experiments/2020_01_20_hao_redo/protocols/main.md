# Protocol for 2020_01_20 Redo of Hao's Hiprfish Staining for Mouse Gut


## materials
- 20" EvaGreen (Biotium; 31000)
- 2" Phusion hot start polymerase master mix (New England Biolabs;
M0536S)
- Tris–EDTA (TE) pH 8 buffer (Ambion; AM9849)
- DNA binding buffer (Zymo Research; D4004-1-L)
- DNA wash buffer (Zymo Research; C1016-50)
- Oligo binding buffer (Zymo Research; D4060-1-40)
- 100-μgcapacitysiliconcolumns(Spin-V;ZymoResearch;D4003-2-48)
- RNA binding buffer (Optional; Zymo Research; R1013-2-100)
- RNA prep buffer (Optional; Zymo Research; R1060-2-100)
- RNA wash buffer (Optional; Zymo Research; R1003-3-24)
- Quick HiScribe T7 polymerase kit (New England Biolabs; E2050S)
- RNasin plus (Promega; N2611)
- Maxima H reverse transcriptase (ThermoScientific; EP0751)
- 10 mM mix of dNTPs (New England Biolabs; N0447S)
- 0.5 M EDTA (Ambion; AM9261)
- 1 N NaOH (VWR; JT5635-2)
- Nuclease-free water (Ambion; AM9932)
- 100% Ethanol (KOPTEC; VWR; 89125-186)
- D/RNAaseFree (VWR; 47751-044)
- 1.5-mL LoBind tubes (Eppendorf; 022431021)
- PCR tubes
- Ethanol (95%, ice cold; 70%)
- Isopropanol (optional; see Step 9)
- 3 M sodium acetate
-
- Sectioned sample in paraffin
- Xylene substitute
- 100% EtOH
- 95% EtOH
- 1x PBS
- Oven at 60˚C
- Wash containers: 50ml falcon tubes or slide rack and coplins jars
- tissue slice formalin fixed on a slide
- lysozyme 100mg/ml in Tris-HCL 10mM pH7.5
- biorad frameseal chambers
- PBS 1x
- EtOH 50%
- Hybridization buffer
  - 20x SSC 10ul
  - Denhardt's sol 10ul
  - Dextran sulfate 20ul
  - Ethylene carbonate 20ul
  - 1% SDS 1ul
  - probe Xul
  - RNAse free ultrapure water 38-Xul
- Wash buffer
  - 5M NaCl 2.15ml
  - 1M HCl 1ml
  - 0.5M EDTA 0.5ml
  - Fill to 50ml with milli Q water

## protocol
### Complex oligo pool Prep
1. Prepare the PCR.
    1. In a 1.7-mL Eppendorf tube, mix the follow- ing: 40 μL 20" Eva Green, 2 μL 200 μM forward primer, 4 μL 100 μM reverse primer, 1 μL of 80 ng/μL complex oligopool, 353 μL nuclease-free water, and 400 μL 2" Phusion hot start polymerase master mix.
    2. Aliquot 50 μL volumes into 16 PCR tubes.
1. Amplify the template.
    1. Run the following protocol on a qPCR machine: (i) 98°C for 3 min, (ii) 98°C for 10 s, (iii) 65°C for 10 s, (iv) 72°C for 15 s, (v) measure the fluorescence of each sample, (vi) repeat 17 times, and (vii) hold at 4˚C.
1. Purify the template.
    1. In a 15-mL Falcon tube, mix the following: 800 μL of the PCR reaction generated in Step 3 and 4 mL of DNA binding buffer.
    2. Pull this mixture across a 100-μg capacity column using either a vacuum manifold or a centrifuge.
    3. Wash the column twice with 300 μL DNA wash buffer, spinning the column in a tabletop centrifuge at maximum speed for 30 s each time.
    4. Elute the template by adding 170 μL nuclease-free water to the column, transferring the column to a fresh 1.7-mL Eppendorf tube, and spinning at maximum speed for 30 s.
    5. Set aside 10 μL of this reaction for quality control.
        1. Measure the dsDNA concentration on a qubit reader.
        1. The concentration should be between 10 and 50 ng/ul
1. In vitro transcription.
    1. In a fresh 1.7-mL Eppendorf tube, mix the following: 160 μL of the in vitro template created in Section 4.2, 176 μL of nuclease-free water, 250 μL of the NTP buffer mix provided with the Quick HiScribe T7 polymerase kit, 25 μL of RNasin Plus, and 25 μL T7 polymer- ase (from the same HiScribe kit).
    2. Incubate the reaction in a 37°C incubator for 12–16 h. Often the reaction is complete after 6–8 h, but it is typically convenient to leave this reaction overnight.
    1. Remove 20 μL for quality control and purify.
        1. Mix 20 μL of the in vitro reaction with 30 μL nuclease-free water, 100 μL RNA binding buffer, and 150 μL 100% ethanol.
        2. Pass across a 100-μg capac- ity spin column in a tabletop centrifuge.
        3. Wash this column once with 400 μL RNA prep buffer, and then twice with 200 μL RNA wash buffer, each time with a 30-s spin at the maximum speed of the tabletop centrifuge.
        4. Elute the RNA with 100 μL nuclease-free water.
        5. Measure RNA concentration on a quibit.
        6. The concentration of the in vitro transcription should be between 0.5 and 2 μg/μL.
1. Reverse transcription.
    1. To the unpurified in vitro transcription, add the following and mix well: 200 μL 10 mM dNTP mix, 120 μL 200 μM forward primer, 240 μL 5" Maxima buffer, 24 μL RNasin Plus, and 24 μL Maxima H# reverse transcriptase.
    2. Incubate In Situ RNA Imaging with MERFISH 27 in a 50°C water bath for 1 h. It is important to use a water bath, not an air incubator, to insure that the temperature of the sample rises to 50°C quickly.
1. Alkaline hydrolysis.
    1. Split the above reaction into two 1.7-mL Eppendorf tubes and add the following to each: 300 μL 0.5 M EDTA and 300 μL 1 N NaOH.
    2. Incubate in a 95°C water bath for 15 min.
1. Purification of ssDNA probe.
    1. Combine the two aliquots above into a single 50-mL Falcon tube and add the following: 4.8 mL Oligo bind- ing buffer and 19.2 mL 100% ethanol.
    2. Mix well and split equally between eight 100-μg capacity spin columns.
    3. Pull the sample across the columns with a vacuum manifold or via centrifugation.
    4. Wash the columns once with 750 μL DNA wash buffer.
    5. Elute the columns using 100 μL of nuclease-free water.
    6. Combine eluates and set aside 10 μL for quality control.
        1. Measure ssDNA concentration using qubit.
        2. Calculate the total mass of ssDNA probes.
1. Concentration of probe.
    1. Add sufficient Sodium acetate to the ssDNA eluate to bring the Sodium acetate concentration to 0.3M.
    2. Add 2 volumes of ice cold ethanol and mix well.
    3. Store on ice or at -20˚C for at least 1hr.
    4. Centrifuge at maximum speed for at least 10min at 0˚C.
    5. Remove supernatant. Fill the tube halfway with 70% EtOH and centrifuge at maximum speed for 2min at 4˚C.
    6. Repeat step 5.
    7. Leave the tube open at room temperature until the fluids have evaporated.
    8. Divide the total mass of ssDNA probes by 33 to generate pmol. Divide the total pmol by the number of oligos in the pool to generate pmol per oligo. Resuspend the concentrated probes to ~0.01μmol/L/oligo.


### Deparaffinization
1. Place slide in 60˚C oven for 1hr.
2. Wash slide in the following sequence
    1. Xylene 15min x2
    3. 100% EtOH 5min

### Tissue permeabilization
1. incubate cells/tissue with lysozyme for 1hr at 37˚C
2. Wash in PBS
3. Store in 50% EtOH

### Hiprfish Staining
1. Prep Wash buffers and heat to 48˚C.
1. Encoding hybridization
    1. Make encoding buffer to 1000ul using 390ul of complex oligo pool.
    1. Pipette hybridization buffer until the tissue is completely covered.
    2. Place the slides in a small sealed chamber containing a significant volume hybridization buffer without dextran sulfate or probes.
    2. Incubate at 46˚C for 12hr.
    2. Wash in wash buffer for 15min at 48˚C. Rinse in EtOH.
1. If using branch probes, incubate slides with hybridization buffer with branch probes for 1hr at 46˚C
1. Readout Hybridization
    1. Incubate hybridization buffer with readout probes for 30min at RT.
    2. Wash in wash buffer for 15min at 48˚C. Rinse in EtOH.
1. Pipette 15ul resin onto slide and cover with coverslip.

Utilities
=========

QUIT contains a number of utilities. Note that these are actually compiled in two separate modules - ``CoreProgs`` contains the bare minimum of commands for the QUIT tests to run, while the actual ``Utils`` modules contains a larger number of useful tools for neuro-imaging pipelines. Their documentation is combined here.

* `qi coil_combine`_
* `qi rfprofile`_
* `qi affine`_
* `qi complex`_
* `qi hdr`_
* `qi kfilter`_
* `qi mask`_
* `qi polyfit/qi polyimg`_
* `qi diff`_
* `qi newimage`_
* `qi select`_

qi coil_combine
---------------

The command implements both the COMPOSER and Hammond methods for coil combination. For COMPOSER, a wrapper script that includes registration and resampling of low resolution reference data to the image data can be found in ``qi composer.sh``.

**Example Command Line**

.. code-block:: bash

    qi coil_combine multicoil_data.nii.gz --composer=composer_reference.nii.gz


Both the input multi-coil file and the reference file must be complex valued. Does not read input from ``stdin``. If a COMPOSER reference file is not specifed, then the Hammond coil combination method is used.

**Outputs**

* ``input_combined.nii.gz`` - The combined complex-valued image.

**Important Options**

* ``--composer, -c``

    Use the COMPOSER method. The reference file should be from a short-echo time reference scan, e.g. UTE or ZTE. If

* ``--coils, -C``

    If your input data is a timeseries consisting of multiple volumes, then use this option to specify the number of coils used in the acquisition. Must match the number of volumes in the reference image. Does not currently work with the Hammond method.


* ``--region, -r``

    The reference region for the Hammond method. Default is an 8x8x8 cube in the center of the acquisition volume.

**References**

- `COMPOSER <http://doi.wiley.com/10.1002/mrm.26093>`_
- `Hammond Method <http://linkinghub.elsevier.com/retrieve/pii/S1053811907009998>`_

qi rfprofile
------------

This utility takes a B1+ (transmit field inhomogeneity) map, and reads an excitation slab profile from ``stdin``. The two are multiplied together along the slab direction (assumed to be Z), to produce a relative flip-angle or B1 map.

**Example Command Line**

.. code-block:: bash

    qi rfprofile b1plus_map.nii.gz output_b1_map.nii.gz < input.json

**Example Command Line**

.. code-block:: json

    {
        "rf_pos" : [ -5, 0, 5],
        "rf_vals" : [[0, 1, 0],
                    [0, 2, 0]]
    }

``rf_pos`` specifies the positions that values of the RF slab have been calculated at, which are specified in ``rf_vals``. Note that ``rf_vals`` is an array of arrays - this allows ``qi rfprofile`` to calculate profiles for multiple flip-angles in a single pass. The units for ``rf_pos`` are the same as image spacing in the header (usually mm). ``rf_vals`` is a unitless fraction, relative to the nominal flip-angle.

These values should be generated with a Bloch simulation. Internally, they are used to create a spline to represent the slab profile. This is then interpolated to each voxel's Z position, and the value multiplied by the input B1+ value at that voxel to produce the output.

**Outputs**

* ``output_b1map.nii.gz`` - The relative flip-angle/B1 map

qi affine
--------

This tool applies simple affine transformations to the header data of an image, i.e. rotations or scalings. It was written because of the inconsistent definitions of co-ordinate systems in pre-clinical imaging. Non-primate mammals are usually scanned prone instead of supine, and are quadrupeds instead of bipeds. This means the definitions of superior/inferior and anterior/posterior are different than in clinical scanning. However, several pre-clinical atlases, e.g. Dorr et al, rotate their data so that the clinical conventions apply. It is hence useful as a pre-processing step to adopt the same co-ordinate system. In addition, packages such as SPM or ANTs have several hard-coded assumptions about their input images that are only appropriate for human brains. It can hence be useful to scale up rodent brains by a factor of 10 so that they have roughly human dimensions.

**Example Command Line**

.. code-block:: bash

    qi affine input_image.nii.gz --scale=10.0 --rotX=90

If no output image is specified, the output will be written back to the input filename.

**Common Options**

- ``--scale, -s``

    Multiply the voxel spacing by a constant factor.

- ``--rotX, --rotY, --rotZ``

    Rotate about the specified axis by the specified number of degrees. Note that currently, each rotation can only be specified once and the order will always be X, Y, then Z.

- ``--offX, --offY, --offZ``

    Add the specified offset to the origin.

- ``--center, -c``

    Set the image origin to be the Center of Gravity of the image.

qi complex
---------

Manipulate complex/real/imaginary/magnitude/phase data. Created because I was fed up with how ``fslcomplex`` works.

**Example Command Line**

.. code-block:: bash

    qi complex -m input_magnitude.nii.gz -p input_phase.nii.gz -R output_real.nii.gz -I output_imaginary.nii.gz

Lower case arguments ``--mag, -m, --pha, -p, --real, -r, --imag, -i, --complex, -x`` are inputs (of which it is only valid to specify certain combinations, complex OR magnitude/phase OR real/imaginary).

Upper case arguments ``--MAG, -M, --PHA, -P, --REAL, -R, --IMAG, -I, --COMPLEX, -X`` are outputs, any or all of which can be specified.

An additional input argument, ``--realimag`` is for Bruker "complex" data, which consists of all real volumes followed by all imaginary volumes, instead of a true complex datatype.

The ``--fixge`` argument fixes the lack of an FFT shift in the slab direction on GE data by multiplying alternate slices by -1. ``--negate`` multiplies the entire volume by -1. ``--double`` reads and writes double precision data instead of floats.

qi hdr
-----

Prints the header of input files as seen by ITK to ``stdout``. Can extract single header fields or print the entirety.

**Example Command Line**

.. code-block:: bash

    qi hdr input_file1.nii.gz input_file2.nii.gz --verbose

Multiple files can be queried at the same time. The ``--verbose`` flag will make sure you can tell which is which.

**Important Options**

If any of the following options are specified, then only those fields will be printed instead of the full header. This is useful if you want to use a header field in a script:
* ``--origin, -o``
* ``--spacing, -S`` - The voxel spacing
* ``--size, -s`` - The matrix size
* ``--voxvol, -v`` - The volume of one voxel

Another useful option is ``--meta, -m``. This will let you query specific image meta-data from the header. You must know the exact name of the meta-data field you wish to obtain.

qi kfilter
---------

MR images often required smoothing or filtering. While this is best done during reconstruction, sometimes it is required as a post-processing step. Instead of filtering by performing a convolution in image space, this tool takes the Fourier Transfrom of input volumes, multiplies k-Space by the specified filter, and transforms back.

**Example Command Line**

.. code-block:: bash

    qi kfilter input_file.nii.gz --filter=Gauss,0.5

**Outputs**

- ``input_file_filtered.nii.gz``

**Important Options**

- ``--filter,-f``

    Specify the filter to use. For all filters below the value \(r\) is the fractional distance from k-Space center, i.e. :math:`r = \sqrt(((k_x / s_x)^2 + (k_y / s_y)^2 + (k_z / s_z)^2) / 3)`. Valid filters are:

    - ``Tukey,a,q``

        A Tukey filter with parameters `a` and `q`. Filter value is 1 for :math:`r < (1 - a)` else the value is :math:`\frac{(1+q)+(1-q)\cos(\pi\frac{r - (1 - a)}{a})}{2}`
    
    - ``Hamming,a,b``

        A Hamming filter, parameters `a` and `b`, value is :math:`a - b\cos(\pi(1+r))`
    
    - ``Gauss,w`` or ``Gauss,x,y,z``

        A Gaussian filter with FWHM specified either isotropically or for each direction independantly.

    - ``Blackman`` or ``Blackman,a``

        A Blackman filter, either with the default parameter of :math:`\alpha=0.16` or the specified :math:`\alpha`. Refer to Wikipedia for the relevant equation.
    
    - ``Rectangle,Dim,Width,Inside,Outside``

        A rectangular or top-hat filter along the specified dimension (must be 0, 1 or 2).
    
    If multiple filters are specified, they are concatenated, *unless* the ``--filter_per_volume`` option is specified.

- ``--filter_per_volume``

    For multiple flip-angle data, the difference in contrast between flip-angles can lead to different amounts of ringing. Hence you may wish to filter volumes with more ringing more heavily. If this option is specified, the number of filters on the command line must match the number of volumes in the input file, and they will be processed in order.

- ``--complex_in`` and ``--complex_out``

    Read / write complex data.

qi mask
------

Implements several different masking strategies. For human data, BET, antsBrainExtraction of 3dSkullStrip are likely better ideas. For pre-clinical data, the strategies below can provide a reasonable mask with some tweaking. There are potentially three stages to generating the mask:

1 - Binary thresholding. If lower or upper thresholds are specified, these are used to separate the image into foreground and background. If neither are specified, then Otsu's method is used to automatically estimate a reasonable threshold value.
2 - (Optional) Run the RATs algorithm
3 - (Optional) Hole-filling

**Example Command Line**

.. code-block:: bash

    qi mask input_image.nii.gz --lower=10 --rats=1200 --fillh=1

In this case an intensity value of 10 will be used as the threshold, RATs will be run with a target volume of 1200 mm^3, and then holes with a radius of 1 voxel will be filled.

**Outputs**

- ``input_image_mask.nii.gz``

**Important Options**

- ``--lower,-l``/``--upper,-u``

    Specify lower and/or upper intensity thresholds. Values below/above these values are set to 0, those inside are set 1. If this option is not specified, Otsu's method will be used to generate a threshold value. If no thresholding is desired, specify ``--lower=0``.

- ``--rats, -r``

    Use the RATs algorithm to remove non-brain tissue. The RATs algorithm uses erode & dilate filters of progressively increasing size until the largest connected component falls below a target size. For rats, target values of around 1000 mm^3 are reasonable.

- ``--fillh, -F``

    Fill holes in the mask up to radius N voxels.

**References**

- `RATs algorithm <http://dx.doi.org/10.1016/j.jneumeth.2013.09.021>`_

qi polyfit/qi polyimg
-------------------

These tools work together to fit Nth order polynomials to images. This is typically used for smoothing a B1 field.

``qi polyfit`` will output the polynomial co-efficients and origin to ``stdout``. ``qi polyimg`` can then read these to generate the polyimage image, using a different image as the reference space. In this way the polynomial image can be created without having to use upsampling.

**Example Command Line**

.. code-block:: bash

    qi polyfit noisy_b1_map.nii.gz --mask=brain_mask.nii.gz --order=8 | qi polyimg hires_t1_image.nii.gz hires_smooth_b1_map.nii.gz --order=8

With the above command-line the output of ``qi polyfit`` is piped directly to the output of ``qi polyimg``. You can instead redirect it to a file with ``>`` and read it in separately. The ``--order`` argument must match between the two commands.

**Important Options**

- ``--order, -o``

    The order of the fitted polynomial. Default is 2 (quadratic)

- ``--mask, -m``

    Only fit the data within a mask. This is usually the brain or only white-matter.

- ``--robust`` (``qi polyimg`` only)

    Use Robust Polynomial Fitting with Huber weights. There is a good discussion of this topic in the Matlab help files.

qi ssfp_bands
-------------

There are several different methods for removing SSFP bands in the literature. Most of them rely on acquiring multiple SSFP images with different phase-increments (also called phase-cycling or phase-cycling patterns). Changing the phase-increments moves the bands to a different location, after which the images can be combined to reduce the banding. The different approaches are discussed further below, but the recommended method is the Geometric Solution which requires complex data.

**Example Command Line**

.. code-block:: bash

    qi ssfpbands ssfp.nii.gz --method=G --2pass --magnitude

The SSFP file must be complex-valued to use the Geometric Solution or Complex Average methods. For the other methods magnitude data is sufficient. Phase-increments should be in opposing pairs, e.g. 180 & 0 degrees, 90 & 270 degrees. These should either be ordered in two blocks, e.g. 180, 90, 0, 270, or alternating, e.g. 180, 0, 90, 270.

**Outputs**

The output filename is the input filename with a suffix that will depend on the method selected (see below).

**Important Options**

- ``--method``

    Choose the band removal method. Choices are:

    - ``G`` Geometric solution. Suffix will be ``GSL`` or ``GSM``
    - `X`` Complex Average. Suffix will be ``CS`` (for Complex Solution)
    - ``R`` Root-mean-square. Suffix will be ``RMS``
    - ``M`` Maximum of magnitudes. Suffix will be ``Max``
    - ``N`` Mean of magnitudes. Suffix will be ``MagMean``

- ``--regularise``

    The Geometric Solution requires regularisation in noisy areas. Available methods are:

    - ``M`` Magnitude regularisation as in original paper
    - ``L`` Line regularisation (unpublished)
    - ``N`` None

    The default is ``L``. If ``L`` or ``M`` are selected, then that character will be appended to the suffix.

- ``--2pass, -2``

    Apply the second-pass energy-minimisation filter from the original paper. Can be likened to smoothing the phase data. If selected will append ``2`` to the suffix.

- ``--alt-order``

    Phase-increments alternate, e.g. 180, 0, 90, 270. The default is the opposite (two blocks), e.g. 180, 90, 0, 270.

- ``--ph-incs``

    Number of phase-increments. The default is 4. If you have multiple phase-increments and (for example) multiple flip-angles, ``qi ssfpbands`` can process them all in one pass.

- ``--ph-order``

    The data order is phase-increment varying fastest, flip-angle slowest. The default is the opposite.

**References**

- `Geometric Solution <http://doi.wiley.com/10.1002/mrm.25098>`_

qi diff
------

Calculates the mean square difference between two images and checks if it is below a tolerance value. Used in the QUIT tests to ensure that calculated parameter maps are close to their baseline values.

**Example Command Line**

.. code-block:: bash

    qi diff --baseline=original.nii --input=calculated.nii --noise=0.01

The command returns the dimensionless noise factor on `stdout`, which is read by the test suite. Note, to make useage clearer, unlike most other QUIT commands all input is specified as arguments.

**Important Options**

- ``--baseline``

    The baseline image. Required.

- ``--image``

    The image to compare to the baseline. Required.

- ``--noise``

    The added noise level.

- ``--tolerance``

    The tolerance is relative to the added noise level (i.e. it is a noise amplification factor).

- ``--abs, -a``

    Use absolute difference instead of fractional difference (i.e. do not divide by the baseline image). Useful when images contain genuine zeros (e.g. off resonance maps).

qi newimage
----------

Creates new images filled with specified patterns. Used for generating test data.

**Example Command Line**

.. code-block:: bash

    qi newimage --size 32,32,32 --grad "0 0.5 1.5" output_image.nii.gz

The file specified on the command line is the *output* file.

**Important Options**

- ``--dims, -d``

    The output dimension. Valid values are 3 and 4.

- ``--size, -s``

    Matrix size of the output image.

- ``--fill, -f``

    Set all voxels in the image to the specified value.

- ``--grad, -g "DIM,LOW,HIGH"``

    Fill voxels with a gradient along the specified dimension, starting at the low value at one edge and finishing at the high value on the other. It is recommended to encase ``DIM,LOW,HIGH`` with quotation marks as they must be passed as a single string to be interpreted properly.

- ``--step, -t "DIM,LOW,HIGH,STEPS"``

    Similar to ``--grad``, but instead of a smooth gradient will with a number of discrete steps.

- ``--wrap, -w``

    Wrap output voxels at the specified value. Useful for simulating phase data.

qi select
---------

Selects a set of volumes from a 4D file and writes them to a new 4D file (a reimplemention of fslselectvols).

**Example Command Line**

.. code-block:: bash

    qi select in_file.nii out_file.nii 2,4,6,8

The last argument is a comma-separated list of the volumes you wish to select.

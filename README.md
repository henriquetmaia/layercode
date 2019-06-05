# layercode
LayerCode: Optical Barcodes for 3D Printed Shapes

In this folder:
├── compressed-dataset
├── hardware-modification
└── matlab
        ├── code
        └── sample-input-images
                ├── successful
                └── failed


The compressed-dataset shows a random view selected
from the shapes explored in our virtual database.
The images are heavily compressed to fit the size limit.

Alongside with the matlab code, we provide sample 
full-resolution images to run the algorithm.

Inside 'successful' we showcase some successful decoding attempts
taken from our database. 

In 'failed' we showcase some difficult LayerCode tags which
our decoding algorithm fails to recover.

To run the decoding algorithm, from the matlab/code folder,
please call

"viewAngleCheck('../sample-input-images/successful');"


Afterwards the results of decoding each folder of images
can be found in the statLog.txt file produced in the
corresponding folder. 
A 1 indicates success, 0 a failed decoding.

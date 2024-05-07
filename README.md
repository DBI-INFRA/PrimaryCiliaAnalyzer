**Requirements**: 
Matlab base package + image processing toolbox (tested with Matlab 2022b)

**Usage:**
  - Download the repository and unzip the sample image
  - Run the Matlab script (double click file and press green arrow in Matlab)
  - Refer to detailed instructions in the script file

**Output:**

Images (maximum intensity projections):
   1) Segmented nuclei (if ShowNucImg == 1)
   2) Cilia candidates (thin red), primary cilia (green)<br>
      Nuclei with (magenta) and without (red) primary cilium<br>
      Nuclei centroids (blue), primary cilia bases (yellow)
   3) Nuclei with valid primary cilium (blue)<br>
      Color-coded main measurement centered on cilia bases

Matlab log:
  - Maximum intensity in nuclei and primary cilia marker channels
  - Total number of nuclei and nuclei with valid primary cilium
  - Primary cilium base, body and average intensity (cilia protein)
  - Primary cilium length (pixels) and base to body intensity ratio

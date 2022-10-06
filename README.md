# CAPs-Toolbox
A toolbox for creating and analyzing Co-Activation Patterns from rs-fMRI data.

## VERSION 0.3 ##

The following scripts are (obliged) to be launched in the sequence:
	1. unfold_and_rename.py
	2. Script_onepop_univr.m

At the moment they are supposed to be run on Linux ( .py) and Windows ( .m ) as on the test computer some other scripts such
as the FSL's ones are installed only in the WSL. All files contained here must be keeped in the same folder for work properly.

DESCRIPTION OF THE SCRIPTS

---- unfold_and_rename.py ----

This script aims to create the FramewiseDisplacement of each subject fMRI image (mov_**.txt), extracting all the volumes and unzipping the
resulting images to .nii format (vol####.nii) in order to give the subsequent Matlab script something readable to its inner SPM12 functions.
After doing this the executable writes a subjectdir.txt textfile that will be used as one input to Script_onepop_univr.m; it also save a copy 
of the original image (Copy_***.nii.gz).

The folders of where subjects are to be specified in the Directories.txt file substituting the existing ones, CAREFUL to maintain the \t spacing
and remember to give the WSL (or directly Linux) folder format and path. The Windows data directory can be also replaced but this time
is mandatory to use the Windows format. It's annoying but for now is the best solution; for sure in future all will run from only one environment.

Folders containing the fMRIs must be all in the same form, in the sense that they must share a prefix.
It is free to use whatever prefix one wants, but it has to be specified in the Directories.txt in the PREFIX line. As one can see the default value is sub-.

After its execution you will find the subjectlist.txt file in the directory containing your data.

---- Script_onepop_univr.m ----

This script is divided in 8 steps (6 from the original TBCAPs Toolbox + 2 brand new).

For now no inputs are requested to the user, but one can modify the ROI location based on where it is saved on the computer.
There are few assumptions on the folders, but it is quickly editable since all variables pointing folders are written as global variables.

The steps are:
	1. Loading the data files (prepareData.m script is launched, it generates the DATASET struct from where the data will be taken
		for the following computations)
	2. Specifying the main parameters
	3. Selecting the frames to analyse   
	4. Consensus clustering (if wished to determine the optimum K)
		4.1 this step has to be inspected as K_opt have to be initialized and there is no default value. 
			Best way is to plot the Qual matrix with a barplot and watch at the highest stability value.
	5. Clustering into CAPs
	6. Computing metrics
	7. Save CAPs nifti
	8. Save CAPs for every subject in folder Results/

######

Notes:

All files must be keeped in the same directory. Also in this directory it should be present the folder of TBCAPs_Toolbox ( https://c4science.ch/source/CAP_Toolbox/ ).
For executing Matlab scripts the Matlab workspace must be set in this folder. At the end the scripts folder will be something like this:
	-|Code
	 |
	 -|univr_CAP_ComputeFD.m
	 -|Script_onepop_univr.m
	 -|unfold_and_rename.py
	 -|Directories.txt
	 -|getCAPS_forEachSubject.m
	 -|prepareData.m
	 -|README.txt
	 -|CAP_Toolbox-master
		

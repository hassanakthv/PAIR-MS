# isoms
*Orbitrap IsoR MS data processing library*

In order to measure the isotopic ratio of each elements (CHNO) in a given sample, you need to convert the .raw files to .mzML file using a converter software such as MSConvertGUI. In the routine platform the isotopic ratio will be determined in MS/MS scans with HCD 50. To convert the file from routine platform one should __add HCD and MS 2__ as filters before converting the .raw to .mzML.

To have a nice information about the samples and experiment, you should provide a .csv file (Experiment Design) with the name of these column name in the same order:

_Mandatory fields_

__File__ : Name of the file with .csv at the end

__Sample__ : Name of the sample as you want to see in the result

__Loading__ : Amount of sample loaded into LC-MS in ng (1000ng = 1ug)

__Start__ : Starting time accoring to retention time, in min (0)

__End__ : Ending time according to retention time, in min

__gC__ : 13C/12C ratio, if it is available

__gN__ : 15N/14N ratio, if it is available

__gH__ : 2H/1H ratio, if it is available

__gO__ : 18O/16O ratio, if it is available
_Followed by experiment information - Note that the column header should be presented in the Experiment Design file but can be left out unfilled_
__Description__ : Information about the experiment
__Aim__ : Goal of the experiment
__Date__ : Date of experiment
__Experimenter__ : Responsible person
__Method__ : Method used in the analysis
__No.Sample__ : Number of samples
__Replicates__ : Number of replicate per each samples
__Other__ : Any other information
__Samples__ : Name of the samples
__Controls__ : Name of the controls
__No.Controls__: Number of controls
__LC.Gradient__: Information about the LC separtion
__Duration__ : Analysis duration time
__Instrument__: Specification on the used instruments
__MS1__ : Specific details about MS scans, e.g. scan range, number of microscans...
__MS2__ : Specific details about MS/MS scans, e.g. scan range, number of microscans...

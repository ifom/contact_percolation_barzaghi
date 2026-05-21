// Example of ImageJ (Fiji) macro used to run StarDist for fluorescent nuclei segmentation
// C. Guidolin

setBatchMode(false)
// setBatchMode(true)


MAIN_FOLDER_RAW = "E:/Chiara/MIXED_POPULATIONS/Leonardo_2024/20240202_mixed_populations_RACDN/";

MAIN_FOLDER_SAVE = "E:/Chiara/MIXED_POPULATIONS/Leonardo_2024_Analysis/MixedPop_RACDN_20240202/";


SAMPLE_NAMES = newArray("R100","C20R80","C20R80_RACDN","C35R65","C35R65_RACDN","C50R50","C50R50_RACDN","C65R35","C65R35_RACDN","C80R20","C80R20_RACDN");

SubfolderList = getFileList(MAIN_FOLDER_RAW); // List of subfolders = raw sample names

//// Check correspondence between file names in Subfolderlist and SAMPLE_NAMES
//for (k = 0; k < SubfolderList.length; k++) {
//print(SubfolderList[k]);
//print(SAMPLE_NAMES[k]);
//}
//waitForUser


for (i = 0; i < SubfolderList.length; i++) {

		SUBFOLDERNAME = SubfolderList[i];
		SUBFOLDER_RAW = MAIN_FOLDER_RAW + SUBFOLDERNAME;
		
		print(SUBFOLDERNAME);
		
		MovieList = getFileList(SUBFOLDER_RAW);  // get the list of files in the subfolder = movies corresponding to different FOVs

		FOV = 0;
		
		File.makeDirectory(MAIN_FOLDER_SAVE + SAMPLE_NAMES[i]);
		
		if (i == 0){ // 100% RAB5 - only channel 2
			
			for (j = 0; j < MovieList.length; j++) {
				
				FILE_NAME = MovieList[j];
				
				FILE_NAME = substring(FILE_NAME,0,lengthOf(FILE_NAME)-4); // remove the ".tif"
				
				print(SAMPLE_NAMES[i]);
				print(FILE_NAME);
				
				FOV = FOV + 1;
							
				print("FOV = " + FOV);
				
				File.makeDirectory(MAIN_FOLDER_SAVE + SAMPLE_NAMES[i] + "/FOV_" + FOV + "/");
				
				open(SUBFOLDER_RAW + FILE_NAME + ".tif");
				
				SUBFOLDER_SAVE_R = MAIN_FOLDER_SAVE + SAMPLE_NAMES[i] + "/FOV_" + FOV + "/R/"; // for channel 2
				File.makeDirectory(SUBFOLDER_SAVE_R);
				
				selectWindow(FILE_NAME + ".tif");
				run("Split Channels");
				
				// channel 1 (Phase Contrast)
				selectWindow("C1-" + FILE_NAME + ".tif");
				close();
				
				// channel 3 (CTRL nuclei - EMPTY)
				selectWindow("C3-" + FILE_NAME + ".tif");
				close();
				
				// channel 2 (RAB5 nuclei)
				selectWindow("C2-" + FILE_NAME + ".tif");
				run("Duplicate...", "duplicate range=101-300");
				selectWindow("C2-" + FILE_NAME + ".tif");
				close();
				selectWindow("C2-" + FILE_NAME + "-1.tif");
				makeRectangle(256, 215, 578, 603);
				run("Specify...", "width=800 height=800 x=112 y=112 slice=1 constrain");
				run("Crop");
				run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'C2-" + FILE_NAME + "-1.tif', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'2.0', 'percentileTop':'98.0', 'probThresh':'0.05', 'nmsThresh':'0.3', 'outputType':'Label Image', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
				selectWindow("Label Image");
				run("Image Sequence... ", "select=" + SUBFOLDER_SAVE_R + " dir=" + SUBFOLDER_SAVE_R + " format=TIFF name=Label_ start=101");
				selectWindow("Label Image");
				close();
				selectWindow("C2-" + FILE_NAME + "-1.tif");
				close();
				
				run("Fresh Start");
				run("Collect Garbage");
			}
			
		}
		
		else { // MIXED POP
			
			for (j = 0; j < MovieList.length; j++) {
				
				FILE_NAME = MovieList[j];
				
				FILE_NAME = substring(FILE_NAME,0,lengthOf(FILE_NAME)-4); // remove the ".tif"
				
				print(SAMPLE_NAMES[i]);
				print(FILE_NAME);
				
				FOV = FOV + 1;
							
				print("FOV = " + FOV);
				
				File.makeDirectory(MAIN_FOLDER_SAVE + SAMPLE_NAMES[i] + "/FOV_" + FOV + "/");
				
				open(SUBFOLDER_RAW + FILE_NAME + ".tif");
				
				SUBFOLDER_SAVE_R = MAIN_FOLDER_SAVE + SAMPLE_NAMES[i] + "/FOV_" + FOV + "/R/"; // for channel 2
				File.makeDirectory(SUBFOLDER_SAVE_R);
				SUBFOLDER_SAVE_C = MAIN_FOLDER_SAVE + SAMPLE_NAMES[i] + "/FOV_" + FOV + "/C/"; // for channel 3
				File.makeDirectory(SUBFOLDER_SAVE_C);
				
				selectWindow(FILE_NAME + ".tif");
				run("Split Channels");
				
				// channel 1 (Phase Contrast)
				selectWindow("C1-" + FILE_NAME + ".tif");
				close();
				
				// channel 2 (RAB5 nuclei)
				selectWindow("C2-" + FILE_NAME + ".tif");
				run("Duplicate...", "duplicate range=101-300");
				selectWindow("C2-" + FILE_NAME + ".tif");
				close();
				selectWindow("C2-" + FILE_NAME + "-1.tif");
				makeRectangle(256, 215, 578, 603);
				run("Specify...", "width=800 height=800 x=112 y=112 slice=1 constrain");
				run("Crop");
				run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'C2-" + FILE_NAME + "-1.tif', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'2.0', 'percentileTop':'98.0', 'probThresh':'0.05', 'nmsThresh':'0.3', 'outputType':'Label Image', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
				selectWindow("Label Image");
				run("Image Sequence... ", "select=" + SUBFOLDER_SAVE_R + " dir=" + SUBFOLDER_SAVE_R + " format=TIFF name=Label_ start=101");
				selectWindow("Label Image");
				close();
				selectWindow("C2-" + FILE_NAME + "-1.tif");
				close();
				
				// channel 3 (CTRL nuclei)
				selectWindow("C3-" + FILE_NAME + ".tif");
				run("Duplicate...", "duplicate range=101-300");
				selectWindow("C3-" + FILE_NAME + ".tif");
				close();
				selectWindow("C3-" + FILE_NAME + "-1.tif");
				makeRectangle(256, 215, 578, 603);
				run("Specify...", "width=800 height=800 x=112 y=112 slice=1 constrain");
				run("Crop");
				run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'C3-" + FILE_NAME + "-1.tif', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'2.0', 'percentileTop':'98.0', 'probThresh':'0.05', 'nmsThresh':'0.3', 'outputType':'Label Image', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
				selectWindow("Label Image");
				run("Image Sequence... ", "select=" + SUBFOLDER_SAVE_C + " dir=" + SUBFOLDER_SAVE_C + " format=TIFF name=Label_ start=101");
				selectWindow("Label Image");
				close();
				selectWindow("C3-" + FILE_NAME + "-1.tif");
				close();
	
				run("Fresh Start");
				run("Collect Garbage");
			}
		}
		
}

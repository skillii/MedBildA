#include <cstdlib>
#include <iostream>
#include <string>

#include "Definitions.h"
#include "VolumeHandler.h"
#include "TVL1.h"
#include "LungSegmentation.h"


#define BUFLEN 100

using std::cout;
using std::endl;

bool bSkipDenoise = false;
char denoisedVolumePath[BUFLEN + 1];



bool ParseCommandLineSettings(int argc, char** argv);
void PrintUsage();

int main(int argc, char** argv)
{
#ifdef _WIN32
  const std::string imagePath = "../../data/VESSEL12_05_256.mhd";
#else
  const std::string imagePath = "../data/VESSEL12_05_256.mhd";
#endif


  VolumeHandler volumeHandler;
  bool success;

  FloatImageType::Pointer denoisedImg;


  if(ParseCommandLineSettings(argc, argv) == false)
  {
	  PrintUsage();
	  return EXIT_FAILURE;
  }

  if(!bSkipDenoise)
  {
	  success = volumeHandler.readVolume(imagePath, true);

	  if (!success)
	  {
		std::cout << "Error reading input volume!" << std::endl;
		return EXIT_FAILURE;
	  }

	  // get the image from the volume handler class
	  FloatImageType::Pointer img = volumeHandler.getHostMem();

	  float min_orig_val,max_orig_val;
	  volumeHandler.getOrigMinMaxValues(min_orig_val, max_orig_val);

	  TVL1 myInstance(img, max_orig_val);

	  denoisedImg = myInstance.Denoise();

	  volumeHandler.setVolume(denoisedImg, min_orig_val, max_orig_val);
	  volumeHandler.writeVolume("denoised_volume.mhd");
  }
  else
  {
	  //do not denoise; used precalculated denoised image instead!
	  success = volumeHandler.readVolume(denoisedVolumePath, true);

	  if (!success)
	  {
		std::cout << "Error reading denoised volume!" << std::endl;
		return EXIT_FAILURE;
	  }

	  denoisedImg = volumeHandler.getHostMem();
  }

  FloatImageType::Pointer thresh_img, lung_img;


  //perform the lung segmentation of the denoised image
  LungSegmentation seg(denoisedImg);

  seg.run();

  //write threshold image (intermediate result of lung-segmentation)
  thresh_img = seg.getThresholdImage();

  volumeHandler.setVolume(thresh_img, 0, 255);
  volumeHandler.writeVolumeBinary("threshold_volume.mhd");



  //write extracted lung
  lung_img = seg.getLungImage();

  volumeHandler.setVolume(lung_img, 0, 255);
  volumeHandler.writeVolumeBinary("lung_volume.mhd");

  
  return EXIT_SUCCESS;
}


bool ParseCommandLineSettings(int argc, char** argv)
{
	strcpy(denoisedVolumePath, "denoised_volume.mhd");

	if(argc <= 1)
		return true;
	else if(argc == 3)
	{
		if(strcmp(argv[1], "-skipdenoise") != 0)
			return false;

		strncpy(denoisedVolumePath, argv[2], BUFLEN);
		denoisedVolumePath[BUFLEN] = 0;
		bSkipDenoise = true;

		return true;
	}
	else
		return false;
}

void PrintUsage()
{
	cout << "miaKUSheet1 [-skipdenoise denoisedimage.mhd]" << endl;
	cout << endl;
	cout << "       -skipdenoise denoisedimage.mhd:" << endl;
	cout << "               skips the denoising process and uses" << endl;
	cout << "               the specified denoised image for further calculations" << endl;
}

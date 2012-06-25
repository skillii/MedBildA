#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "Definitions.h"
#include "VolumeHandler.h"
#include "TubeDetection.h"
#include "CenterlineExtraction.h"


#define BUFLEN 100

#define READ_RESULT 1

using std::cout;
using std::endl;



bool ParseCommandLineSettings(int argc, char** argv);
void PrintUsage();

int main(int argc, char** argv)
{
#ifdef _WIN32
    const std::string imagePath = "../../data/VESSEL12_05_256.mhd";
    const std::string lungSegmentPath = "../../data/lung_volume.mhd";
#else
    const std::string imagePath = "../data/VESSEL12_05_256.mhd";
    const std::string lungSegmentPath = "../data/lung_volume.mhd";
#endif


    VolumeHandler volumeHandler;
    bool success;


    if(ParseCommandLineSettings(argc, argv) == false)
    {
      PrintUsage();
      return EXIT_FAILURE;
    }

    success = volumeHandler.readVolume(imagePath, true);

    if (!success)
    {
    std::cout << "Error reading input volume!" << std::endl;
    return EXIT_FAILURE;
    }

    // get the image from the volume handler class
    FloatImageType::Pointer raw = volumeHandler.getHostMem();

    float min_orig_val,max_orig_val;
    volumeHandler.getOrigMinMaxValues(min_orig_val, max_orig_val);

    std::cout << "min orig value: " << min_orig_val << std::endl;
    std::cout << "max orig value: " << max_orig_val << std::endl;

    success = volumeHandler.readVolume(lungSegmentPath, true);

    if (!success)
    {
    std::cout << "Error reading input volume!" << std::endl;
    return EXIT_FAILURE;
    }

    // get the image from the volume handler class
    FloatImageType::Pointer lungSegment = volumeHandler.getHostMem();


    if(READ_RESULT == 0)
    {
      TubeDetection tubeDetector(lungSegment, raw);

      tubeDetector.cropLung();
      tubeDetector.buildImagePyramid();
      tubeDetector.calcHessian();
      tubeDetector.calcGradients();

      tubeDetector.calcMedialness();
      tubeDetector.calcMaxMedialness();

      unsigned i,j;
      char path[101];

      std::cout << "Writing files.." << std::endl;

      ///////////////////////   save image pyramid
      std::vector<FloatImageType::Pointer> imgPyramid = tubeDetector.getImgPyramid();


      for(i = 0; i < imgPyramid.size(); i++)
      {
        sprintf(path, "pyramid/%d.mhd", i);

        volumeHandler.setVolume(imgPyramid.at(i), min_orig_val, max_orig_val);
        volumeHandler.writeVolume(path);
      }

    ///////////////////////   save gradient images

      for (j = 1; j <= 3; j++)
      {
        std::vector<FloatImageType::Pointer> grad = tubeDetector.getGradientImage(j-1);


        for(i = 0; i < grad.size(); i++)
        {
          sprintf(path, "gradients/%d/%d.mhd",j, i);

          volumeHandler.setVolume(grad.at(i), min_orig_val, max_orig_val);
          volumeHandler.writeVolume(path);
        }
      }


    /////////////////////// save Eigenvector images
      for (j = 0; j < 3; j++)
      {
        FloatImageType::Pointer ev = tubeDetector.getEVImage(j);


        sprintf(path, "EV/%d.mhd",j);

        volumeHandler.setVolume(ev, min_orig_val, max_orig_val);
        volumeHandler.writeVolume(path);
      }



    ///////////////////////   save medialness pyramid
      std::vector<FloatImageType::Pointer> medialness = tubeDetector.getMedialnessImages();

      float lo,hi;

      while(1) {
      std::cout << "enter min value: " << std::endl;
      std::cin >> lo;

      std::cout << "enter max value: " << std::endl;
      std::cin >> hi;

      for(i = 0; i < medialness.size(); i++)
      {
          sprintf(path, "medialness/%d.mhd", i);

          volumeHandler.setVolume(medialness.at(i), lo, hi);
           volumeHandler.writeVolume(path);
      }

       //////////////////////// save max medialness image

         sprintf(path, "medialness/max.mhd");
         volumeHandler.setVolume(tubeDetector.getMaxMedialnessImage(),min_orig_val, max_orig_val);
         volumeHandler.writeVolume(path);
      }
    }

    CenterlineExtraction ce;


    int result = ce.readSavedFiles();

    if(result == EXIT_SUCCESS)
    {
        ce.performNonMaximaSurpression();
        ce.performReconnection();

        FloatImageType::Pointer medialness = ce.getMedialnessImage();

        char path[101];

        sprintf(path, "medialness/nonmaxima.mhd");

        volumeHandler.setVolume(medialness, 0, 100);
        volumeHandler.writeVolume(path);


        FloatImageType::Pointer centerline = ce.getCenterlineImage();


        sprintf(path, "medialness/centerline.mhd");

        volumeHandler.setVolume(centerline, 0, 100);
        volumeHandler.writeVolume(path);


    }





    return EXIT_SUCCESS;
}


bool ParseCommandLineSettings(int argc, char** argv)
{
  return true;
}

void PrintUsage()
{
	cout << "miaKUSheet1 [-skipdenoise denoisedimage.mhd]" << endl;
	cout << endl;
	cout << "       -skipdenoise denoisedimage.mhd:" << endl;
	cout << "               skips the denoising process and uses" << endl;
	cout << "               the specified denoised image for further calculations" << endl;
}

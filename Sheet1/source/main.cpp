#include <cstdlib>
#include <iostream>

#include "Definitions.h"
#include "VolumeHandler.h"

int main(int argc, char** argv)
{
  std::cout << "Hello World!" << std::endl;

#ifdef _WIN32
  const std::string imagePath = "../../data/VESSEL12_05_256.mhd";
#else
  const std::string imagePath = "../data/VESSEL12_05_256.mhd";
#endif
 
  VolumeHandler volumeHandler;
  bool success = volumeHandler.readVolume(imagePath, true);

  if (!success)
  {
    std::cout << "Error reading input volume!" << std::endl;
    return EXIT_FAILURE;
  }

  // get the image from the volume handler class
  FloatImageType::Pointer img = volumeHandler.getHostMem();
	
  float min_orig_val,max_orig_val;
  volumeHandler.getOrigMinMaxValues(min_orig_val, max_orig_val);

  std::cout << "min_ and max_orig_val: " << min_orig_val << " ; " << max_orig_val << std::endl;

  // invert picture
  itk::ImageRegionIterator<FloatImageType> it(img, img->GetLargestPossibleRegion().GetSize());
  
  while(!it.IsAtEnd())
  {
    it.Set(1.0f - it.Get());
    ++it;
  }

  volumeHandler.setVolume(img, min_orig_val, max_orig_val);
  volumeHandler.writeVolumeBinary("inverted_volume.mhd");
  
  return EXIT_SUCCESS;
}

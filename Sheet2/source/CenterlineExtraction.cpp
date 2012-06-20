/*
 * CenterlineExtraction.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: tom
 */


#include "CenterlineExtraction.h"
#include "TubeDetection.h"
#include "Numericshelper.h"
#include <cmath>


const float CenterlineExtraction::threshold_low = 0.1f;
const float CenterlineExtraction::threshold_high = 0.4f;

int CenterlineExtraction::readSavedFiles()
{
    bool success = volumeHandler.readVolume("EV/0.mhd", true);

    if (!success)
    {
      std::cout << "Error reading input volume EV/0.mhd!" << std::endl;
      return EXIT_FAILURE;
    }

    // get the image from the volume handler class
    this->eigenvector[0] = volumeHandler.getHostMem();

    volumeHandler.getOrigMinMaxValues(min_vals_eigenvector[0], max_vals_eigenvector[0]);

    success = volumeHandler.readVolume("EV/1.mhd", true);

    if (!success)
    {
      std::cout << "Error reading input volume EV/1.mhd!" << std::endl;
      return EXIT_FAILURE;
    }

    // get the image from the volume handler class
    this->eigenvector[1] = volumeHandler.getHostMem();

    volumeHandler.getOrigMinMaxValues(min_vals_eigenvector[1], max_vals_eigenvector[1]);

    success = volumeHandler.readVolume("EV/2.mhd", true);

    if (!success)
    {
      std::cout << "Error reading input volume EV/2.mhd!" << std::endl;
      return EXIT_FAILURE;
    }

    // get the image from the volume handler class
    this->eigenvector[2] = volumeHandler.getHostMem();

    volumeHandler.getOrigMinMaxValues(min_vals_eigenvector[2], max_vals_eigenvector[2]);


    success = volumeHandler.readVolume("medialness/max.mhd", true);

    if (!success)
    {
      std::cout << "Error reading input volume medialness/max.mhd!" << std::endl;
      return EXIT_FAILURE;
    }

    // get the image from the volume handler class
    this->maxMedialnessOverScales = volumeHandler.getHostMem();

    volumeHandler.getOrigMinMaxValues(this->min_val_maxMedialness, this->max_val_maxMedialness);

    return EXIT_SUCCESS;
}


void printVector(itk::Vector<float, 3> v)
{
    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
}

void CenterlineExtraction::performNonMaximaSurpression()
{
    unsigned x,y,z;
    int i;

    unsigned int x_size;
    unsigned int y_size;
    unsigned int z_size;

    FloatImageType::IndexType index;

    std::cout << "Starting non-Maxima supression" << std::endl;

    x_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(0);
    y_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(1);
    z_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(2);

    float current_medialness_value;
    float neighbor_medialness_value;

    itk::Vector<float, 3> v;
    itk::Vector<float, 3> v_i[2];

    float alpha_step = 2*M_PI/TubeDetection::alpha_steps;
    float alpha = 0;

    for(x = 0; x < x_size; x++)
    {
      index[0] = x;

      if(x % 10 == 0)
        std::cout << "x " << x << std::endl;

      for(y = 0; y < y_size; y++)
      {
        index[1] = y;

        for(z = 0; z < z_size; z++)
        {
          index[2] = z;


          current_medialness_value = maxMedialnessOverScales->GetPixel(index);

          if(current_medialness_value < threshold_low)
          {
              maxMedialnessOverScales->SetPixel(index,0);
              continue;
          }


          //Calculate other 2 eigenvector here and perform non-maxima surpression
          v[0] = this->eigenvector[0]->GetPixel(index);
          v[1] = this->eigenvector[1]->GetPixel(index);
          v[2] = this->eigenvector[2]->GetPixel(index);

          v.Normalize();

          if(v[2] < 0.5)
          {
            v_i[0][0] = -v[1];
            v_i[0][1] = v[0];
            v_i[0][2] =  0;

            v_i[0].Normalize();
          }
          else
          {
            v_i[0][0] = -v[2];
            v_i[0][1] = 0;
            v_i[0][2] = v[0];

            v_i[0].Normalize();
          }


          v_i[1] = itk::CrossProduct(v,v_i[0]);
          v_i[1].Normalize();


          //Check out the neighborinos
          for(i = 0; i < TubeDetection::alpha_steps; i++, alpha += alpha_step)
          {
             v = v_i[0]*cosf(alpha) + v_i[1]*sinf(alpha);

             neighbor_medialness_value = NumericsHelper::trilinearInterp(maxMedialnessOverScales,
                                                                         v[0] + x, v[1] + y, v[2] + z);


             if(current_medialness_value < neighbor_medialness_value)
             {
               //Not a local maxima! Set value at index to 0
               maxMedialnessOverScales->SetPixel(index,0);
               break;
             }
          }

        }
      }
    }


}



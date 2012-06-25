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


const float CenterlineExtraction::threshold_low = 0.04f;
const float CenterlineExtraction::threshold_high = 0.15f;
const int CenterlineExtraction::skippixels = 3;

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


void CenterlineExtraction::performReconnection()
{
    FloatImageType::IndexType index;

    itk::Vector<float, 3> v;


    std::cout << "starting reconnection: " << std::endl;

    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(maxMedialnessOverScales);
    duplicator->Update();
    centerlineImage = duplicator->GetOutput();


    SetAllPixels(centerlineImage, 0.0f);



    while(!reconnectionQueue.empty())
    {
        index = reconnectionQueue.front();
        reconnectionQueue.pop();

        v[0] = this->eigenvector[0]->GetPixel(index);
        v[1] = this->eigenvector[1]->GetPixel(index);
        v[2] = this->eigenvector[2]->GetPixel(index);

        //printVector(v);

        performReconnection(index, v);
        //performReconnection(index, -v);
    }

    std::cout << "finished reconnection: " << std::endl;
}

void CenterlineExtraction::performReconnection(FloatImageType::IndexType index, itk::Vector<float, 3> direction)
{
    int x,y,z;
    unsigned int x_size;
    unsigned int y_size;
    unsigned int z_size;

    float pixelvalue;
    float initialvalue;

    int pathlength = 0;
    static int total_pathlength = 0;
    static int count = 0;

    direction.Normalize();


    itk::Vector<float, 3> pos;


    x_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(0);
    y_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(1);
    z_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(2);


    initialvalue = maxMedialnessOverScales->GetPixel(index);
    //initialvalue = 1.0f;

    centerlineImage->SetPixel(index, initialvalue);

    //centerlineImage->SetPixel(index, initialvalue);

    pos[0] = index[0];
    pos[1] = index[1];
    pos[2] = index[2];

    int skipped_positions = 0;
    std::queue<FloatImageType::IndexType> skipped_queue;
    itk::Vector<float, 3> old_direction = direction;

    FloatImageType::IndexType nearestNeighbour = index;

    while(1)
    {
        direction[0] = this->eigenvector[0]->GetPixel(nearestNeighbour);
        direction[1] = this->eigenvector[1]->GetPixel(nearestNeighbour);
        direction[2] = this->eigenvector[2]->GetPixel(nearestNeighbour);

        direction.Normalize();


        if(abs(direction * old_direction) > cos(90))
        {
         if(direction * old_direction < 0)
           direction = -direction;
        }
        else
          direction = old_direction;


        pos = pos + direction;
        old_direction = direction;



        x = static_cast<unsigned>(floor(pos[0]));
        y = static_cast<unsigned>(floor(pos[1]));
        z = static_cast<unsigned>(floor(pos[2]));

        //boundary check:
        if(x < 0 || x >= (int)x_size || y < 0 || y >= (int)y_size || z < 0 || z >= (int)z_size)
            break;


        nearestNeighbour[0] = x; nearestNeighbour[1] = y; nearestNeighbour[2] = z;


        pixelvalue = maxMedialnessOverScales->GetPixel(nearestNeighbour);
        //pixelvalue = NumericsHelper::trilinearInterp(maxMedialnessOverScales, pos[0], pos[1], pos[2]);

        if(pixelvalue > threshold_low)
        {
            if(initialvalue > centerlineImage->GetPixel(nearestNeighbour))
                centerlineImage->SetPixel(nearestNeighbour, initialvalue);

            pathlength++;

            while(skipped_queue.size() != 0)
            {
              std::cout << "skipping(pathlength=" << pathlength++ << ")" << std::endl;
              centerlineImage->SetPixel(skipped_queue.front(), initialvalue);
              skipped_queue.pop();
            }

            skipped_positions = 0;
        }
        else
        {
            /*itk::Vector<float, 3> neighbourdir;

            if(getMaxNeighbour(nearestNeighbour, neighbourdir) > threshold_low)
            {
                if(neighbourdir*direction > -0.5)   //don't go back!
                {
                    centerlineImage->SetPixel(nearestNeighbour, 1);
                    pos = pos + neighbourdir;
                    direction = neighbourdir;
                    std::cout << "----------------------------" << std::endl;
                    std::cout << "----------------------------" << std::endl;
                    std::cout << "----------------------------" << std::endl;
                    std::cout << "----------------------------" << std::endl;
                    pathlength++;
                    continue;
                }
            }*/

            skipped_positions++;

            if(skipped_positions > skippixels)
                break;
            else
                skipped_queue.push(nearestNeighbour);
        }
    }


    total_pathlength += pathlength;
    count++;
    std::cout << "total length of reconnection: " << pathlength << std::endl;
    std::cout << "mean pathlength: " << static_cast<float>(total_pathlength) / count << std::endl;
}

float CenterlineExtraction::getMaxNeighbour(FloatImageType::IndexType index, itk::Vector<float, 3> &direction)
{
    float value;

    int x,y,z;

    FloatImageType::IndexType ind;

    FloatImageType::IndexType maxInd;

    float max = 0;

    for(x = -1; x <= 1; x++)
    {
        for(y = -1; y <= 1; y++)
        {
            for(z = -1; z <= 1; z++)
            {
                if(x == 0 && y == 0 && z == 0)
                    continue;

                ind[0] = index[0] + x;
                ind[1] = index[1] + y;
                ind[2] = index[2] + z;

                value = maxMedialnessOverScales->GetPixel(ind);

                if(value >= max)
                {
                    max = value;
                    maxInd = ind;
                }
            }
        }
    }

    direction[0] = maxInd[0] - ind[0];
    direction[1] = maxInd[1] - ind[1];
    direction[2] = maxInd[2] - ind[2];

    return max;
}



void CenterlineExtraction::SetAllPixels(FloatImageType::Pointer img, float value)
{
    unsigned int x, y, z;

    unsigned int x_size;
    unsigned int y_size;
    unsigned int z_size;

    FloatImageType::IndexType index;


    x_size = img->GetLargestPossibleRegion().GetSize(0);
    y_size = img->GetLargestPossibleRegion().GetSize(1);
    z_size = img->GetLargestPossibleRegion().GetSize(2);


    for(x = 0; x < x_size; x++)
    {
        index[0] = x;

        for(y = 0; y < y_size; y++)
        {
            index[1] = y;

            for(z = 0; z < z_size; z++)
            {
                index[2] = z;

                img->SetPixel(index, value);
            }
        }
    }
}

void CenterlineExtraction::performNonMaximaSurpression()
{
    unsigned x,y,z;
    int i;

    int isMax;

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

          isMax = 1;

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
               isMax = 0;
               break;
             }
          }

          if(isMax && current_medialness_value > threshold_high)
              reconnectionQueue.push(index);

        }
      }
    }


}



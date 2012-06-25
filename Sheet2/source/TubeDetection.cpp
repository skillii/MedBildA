/*
 * TubeDetection.cpp
 *
 *  Created on: Jun 16, 2012
 *      Author: tom
 */

#include <cmath>
#include <iostream>
#include "TubeDetection.h"
#include "VolumeHandler.h"
#include "Numericshelper.h"



#define START_LEVEL (0) //number of scale-levels to skip

const int TubeDetection::scale_levels = SCALE_LEVELS;
const float TubeDetection::tube_r[SCALE_LEVELS] = {1, 1, 1};
const int TubeDetection::alpha_steps = 8;



TubeDetection::TubeDetection(FloatImageType::Pointer lungSegment, FloatImageType::Pointer raw) {
    this->lungSegment = lungSegment;
    this->raw = raw;
    this->lung = raw;
}

TubeDetection::~TubeDetection() {
}


void TubeDetection::cropLung()
{
    unsigned x,y,z;
    FloatImageType::IndexType index;
    unsigned int x_size = raw->GetLargestPossibleRegion().GetSize(0);
    unsigned int y_size = raw->GetLargestPossibleRegion().GetSize(1);
    unsigned int z_size = raw->GetLargestPossibleRegion().GetSize(2);

    std::cout << "Starting crop lung" << std::endl;

    for(x = 0; x < x_size ; x++)
    {
        index[0]    = x;


        for(y = 0; y < y_size; y++)
        {
            index[1] = y;


            for(z = 0; z < z_size; z++)
            {
                //Indices for gradient, forward differences
                index[2] = z;


                lung->SetPixel(index, raw->GetPixel(index)*lungSegment->GetPixel(index));
                //lung->SetPixel(index, raw->GetPixel(index));
            }
        }
    }

    std::cout << "Crop lung done" << std::endl;
}


void TubeDetection::buildImagePyramid()
{
    PyramidFilter::Pointer p_filter = PyramidFilter::New();

    std::cout << "Starting pyramid construction" << std::endl;

    p_filter->SetInput(lung);
    p_filter->SetStartingShrinkFactors(2);
    p_filter->SetNumberOfLevels(scale_levels);

    p_filter->Update();


    int levels = p_filter->GetNumberOfLevels();

    for(int i = 0; i < levels; i++)
    {
        if(i >= (levels-START_LEVEL))
            continue;
        imagePyramid.push_back(p_filter->GetOutput(i));
    }

    std::cout << "Finished pyramid construction" << std::endl;
}


void TubeDetection::calcHessian()
{
    std::cout << "Starting hessian calculation" << std::endl;



    for(unsigned i = 0; i < this->imagePyramid.size(); i++)
    {
        HessianFilter::Pointer h_filter = HessianFilter::New();

        h_filter->SetInput(imagePyramid.at(i));

        h_filter->SetSigma(1);

        h_filter->Update();

        hessianImages.push_back(h_filter->GetOutput());

    }


    std::cout << "Finished hessian calculation" << std::endl;
}

void TubeDetection::calcGradients()
{
    std::cout << "Starting gradient calculation" << std::endl;

    for(unsigned i = 0; i < this->imagePyramid.size(); i++)
    {
        std::cout << "Calculating gradient for pyramid level " << i << std::endl;

        GradientImageFilter::Pointer g_filter = GradientImageFilter::New();

        g_filter->SetInput(imagePyramid.at(i));

        g_filter->Update();


        GradientImageFilter::OutputImagePointer grad = GradientImageFilter::OutputImageType::New();

        grad = g_filter->GetOutput();


        //convert gradient to float image types

        FloatImageType::Pointer gradX, gradY, gradZ;

        //allocate new images for the gradients

        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage(imagePyramid.at(i));
        duplicator->Update();

        gradX = duplicator->GetOutput();

        duplicator = DuplicatorType::New();
        duplicator->SetInputImage(imagePyramid.at(i));
        duplicator->Update();

        gradY = duplicator->GetOutput();

        duplicator = DuplicatorType::New();
        duplicator->SetInputImage(imagePyramid.at(i));
        duplicator->Update();

        gradZ = duplicator->GetOutput();


        GradientImageFilter::OutputImageType::IndexType index;
        FloatImageType::IndexType index_f;

        unsigned int x_size;
        unsigned int y_size;
        unsigned int z_size;


        x_size = grad->GetLargestPossibleRegion().GetSize(0);
        y_size = grad->GetLargestPossibleRegion().GetSize(1);
        z_size = grad->GetLargestPossibleRegion().GetSize(2);

        unsigned int x,y,z;
        for(x = 0; x < x_size ; x++)
        {
            index[0] = x;
            index_f[0] = x;

            for(y = 0; y < y_size; y++)
            {
                index[1] = y;
                index_f[1] = y;

                for(z = 0; z < z_size; z++)
                {
                    index[2] = z;
                    index_f[2] = z;

                    GradientImageFilter::OutputImagePixelType gradient = grad->GetPixel(index);


                    gradX->SetPixel(index_f, gradient[0]);
                    gradY->SetPixel(index_f, gradient[1]);
                    gradZ->SetPixel(index_f, gradient[2]);

                }
            }
        }
        gradientX.push_back(gradX);
        gradientY.push_back(gradY);
        gradientZ.push_back(gradZ);
    }

    std::cout << "Finished gradient calculation" << std::endl;
}

void TubeDetection::allocateEigenvectorImage()
{
    FloatImageType::IndexType start;
    FloatImageType::SizeType size;
    FloatImageType::RegionType region;

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    size[0] = imagePyramid.back()->GetLargestPossibleRegion().GetSize(0);
    size[1] = imagePyramid.back()->GetLargestPossibleRegion().GetSize(1);
    size[2] = imagePyramid.back()->GetLargestPossibleRegion().GetSize(2);

    eigenvector[0] = FloatImageType::New();
    eigenvector[1] = FloatImageType::New();
    eigenvector[2] = FloatImageType::New();


    region.SetSize(size);
    region.SetIndex(start);



    eigenvector[0]->SetRegions(region);
    eigenvector[1]->SetRegions(region);
    eigenvector[2]->SetRegions(region);


    eigenvector[0]->Allocate();
    eigenvector[1]->Allocate();
    eigenvector[2]->Allocate();
}

void TubeDetection::eigenValueDecomposition(int scale, HessianFilter::OutputImageType::IndexType index, float V[3][3], float d[3])
{
    HessianFilter::OutputImagePixelType current_hessian;
    float *float_ptr;
    double *double_ptr;
    float a[3][3];
    double* a_double;

    current_hessian = hessianImages.at(scale)->GetPixel(index);


    //current_hessian.FixedArray(a_double);
    a_double = current_hessian.GetDataPointer();

    int u=0;
    for(float_ptr = &a[0][0], double_ptr = a_double, u=0; u<6; u++, float_ptr++, double_ptr++)
        *float_ptr = (float)*double_ptr;


    float V_unsorted[3][3];
    float d_unsorted[3];

    NumericsHelper::calculateEigenDecompositionSymmetric3x3(a,V_unsorted,d_unsorted);



    //sort Eigenvalues: first one is the smallest one:

    int i,j;
    int min_ev;
    NumericsHelper::min(d_unsorted, 3, &min_ev);

    d[0] = d_unsorted[min_ev];
    V[0][0] = V_unsorted[min_ev][0];
    V[0][1] = V_unsorted[min_ev][1];
    V[0][2] = V_unsorted[min_ev][2];



    for(i = 0,j = 1; i < 3; i++)
    {
        if(i != min_ev)
        {
            V[j][0] = V_unsorted[i][0];
            V[j][1] = V_unsorted[i][1];
            V[j][2] = V_unsorted[i][2];
            d[j] = d_unsorted[i];
            j++;
        }
    }
}

void TubeDetection::calcMedialness()
{
    unsigned x,y,z;
    HessianFilter::OutputImageType::IndexType index;
    unsigned int x_size;
    unsigned int y_size;
    unsigned int z_size;



    float V[3][3];
    float d[3];

    FloatImageType::Pointer medialnessimage;

    std::cout << "Starting calculating medialness for level: ";


    for(unsigned i = 0; i < this->imagePyramid.size(); i++)
    {
        std::cout << i << " " << std::endl;

        //medialnessimage = imagePyramid.at(i)->Clone();
        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage(imagePyramid.at(i));
        duplicator->Update();

        medialnessimage = duplicator->GetOutput();

        x_size = imagePyramid.at(i)->GetLargestPossibleRegion().GetSize(0);
        y_size = imagePyramid.at(i)->GetLargestPossibleRegion().GetSize(1);
        z_size = imagePyramid.at(i)->GetLargestPossibleRegion().GetSize(2);

        for(x = 0; x < x_size ; x++)
        {
            index[0] = x;

            if(x % 10 == 0)
              std::cout << "x = " << x << std::endl;

            for(y = 0; y < y_size; y++)
            {
                index[1] = y;

                for(z = 0; z < z_size; z++)
                {
                    //Indices for gradient, forward differences
                    index[2] = z;


                    eigenValueDecomposition(i, index, V, d);

                    float medialness = calcMedialness(i, x, y, z, d, V);

                    medialnessimage->SetPixel(index, medialness);
                }
            }
        }

        medialnessImages.push_back(medialnessimage);
    }

    std::cout << std::endl << "Finished calculation of medialness" << std::endl;
}

void TubeDetection::calcMaxMedialness()
{
    unsigned x,y,z;

    unsigned int x_size;
    unsigned int y_size;
    unsigned int z_size;

    float V[3][3];
    float d[3];

    FloatImageType::Pointer max_medialnessimage;

    FloatImageType::IndexType index;

    float scale_index[3];
    HessianFilter::OutputImageType::IndexType indexi;

    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(medialnessImages.at(medialnessImages.size() - 1));
    duplicator->Update();

    max_medialnessimage = duplicator->GetOutput();

    allocateEigenvectorImage();


    std::cout << "Starting calculating maximum of medialness." << std::endl;

    x_size = medialnessImages.at(medialnessImages.size() - 1)->GetLargestPossibleRegion().GetSize(0);
    y_size = medialnessImages.at(medialnessImages.size() - 1)->GetLargestPossibleRegion().GetSize(1);
    z_size = medialnessImages.at(medialnessImages.size() - 1)->GetLargestPossibleRegion().GetSize(2);

    float current_medialness_value;
    float max_medialness_value;
    int  max_medialness_scale_level;

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

          //This is the value we obtain from the largest resolution
          max_medialness_value = max_medialnessimage->GetPixel(index);
          max_medialness_scale_level = medialnessImages.size() - 1;

          //Now we go through all scale levels and look for a maximum..
          for(unsigned i = 0; i < this->medialnessImages.size() - 1; i++)
          {
            scale_index[0] = static_cast<float>(x) / powf(2,medialnessImages.size() - 1 - i);
            scale_index[1] = static_cast<float>(y) / powf(2,medialnessImages.size() - 1 - i);
            scale_index[2] = static_cast<float>(z) / powf(2,medialnessImages.size() - 1 - i);

            current_medialness_value = NumericsHelper::trilinearInterp(medialnessImages.at(i), scale_index[0], scale_index[1], scale_index[2]);

            if(current_medialness_value > max_medialness_value)
            {
              max_medialness_value = current_medialness_value;
              max_medialness_scale_level = i;
            }

          }

          //Now write the maximum into the image again
          max_medialnessimage->SetPixel(index, max_medialness_value);



          //now calculate and save the EV of that pixel:
          indexi[0] = static_cast<int>(roundf(static_cast<float>(x) / powf(2,medialnessImages.size() - 1 - max_medialness_scale_level)));
          indexi[1] = static_cast<int>(roundf(static_cast<float>(y) / powf(2,medialnessImages.size() - 1 - max_medialness_scale_level)));
          indexi[2] = static_cast<int>(roundf(static_cast<float>(z) / powf(2,medialnessImages.size() - 1 - max_medialness_scale_level)));


          eigenValueDecomposition(max_medialness_scale_level, indexi, V, d);


          eigenvector[0]->SetPixel(index, V[0][0]);
          eigenvector[1]->SetPixel(index, V[0][1]);
          eigenvector[2]->SetPixel(index, V[0][2]);

        }
      }
    }

   this->maxMedialnessOverScales = max_medialnessimage;
   std::cout << "Finished calculating maximum of medialness." << std::endl;
}

float TubeDetection::calcMedialness(unsigned level, unsigned x, unsigned y, unsigned z, float ew[3], float ev[3][3])
{
    int i;


    float alpha_step = 2*M_PI/alpha_steps;
    float alpha = 0;

    float c[alpha_steps];

    FloatImageType::Pointer img = imagePyramid.at(level);

    itk::Vector<float, 3> v, grad;
    itk::Vector<float, 3> v_i[2];


    v_i[0].SetElement(0, ev[1][0]);
    v_i[0].SetElement(1, ev[1][1]);
    v_i[0].SetElement(2, ev[1][2]);

    v_i[1].SetElement(0, ev[2][0]);
    v_i[1].SetElement(1, ev[2][1]);
    v_i[1].SetElement(2, ev[2][2]);



    float c_sum = 0, c_var_sum = 0, c_mean, c_var;


    for(i = 0; i < alpha_steps; i++, alpha += alpha_step)
    {
        v = v_i[0]*cosf(alpha) + v_i[1]*sinf(alpha);


        float g_x,g_y,g_z;

        g_x = v[0]*tube_r[level] + x;
        g_y = v[1]*tube_r[level] + y;
        g_z = v[2]*tube_r[level] + z;



        grad[0] = NumericsHelper::trilinearInterp(gradientX.at(level), g_x, g_y, g_z);
        grad[1] = NumericsHelper::trilinearInterp(gradientY.at(level), g_x, g_y, g_z);
        grad[2] = NumericsHelper::trilinearInterp(gradientZ.at(level), g_x, g_y, g_z);

        c[i] = -grad*v;

        if(c[i] < 0)
            c[i] = 0;

        c_sum += c[i];
    }

    c_mean = c_sum / alpha_steps;
    //calculate variance:

    for(i = 0; i < alpha_steps; i++, alpha += alpha_step)
    {
        c_var_sum += powf((c[i] - c_mean),2);
    }

    c_var = c_var_sum / alpha_steps;




    //calc symmetry based medialness:

    float s;

    if(c_mean > 0)
        s = 1 - c_var/powf(c_mean,2);
    else
        s = 0;

    s = s*c_mean;


    FloatImageType::IndexType index;

    index[0] = x;
    index[1] = y;
    index[2] = z;

    grad[0] = gradientX.at(level)->GetPixel(index);
    grad[1] = gradientY.at(level)->GetPixel(index);
    grad[2] = gradientZ.at(level)->GetPixel(index);

    float rc = grad.GetNorm();



    float final = s-0.5*rc;

    if(final < 0)
        final = 0;

    //return final; DUMMY!!!

    return final;
}


/*void TubeDetection::calcNonMaximumSupression()
{
    unsigned x,y,z;

    unsigned int x_size;
    unsigned int y_size;
    unsigned int z_size;


    FloatImageType::IndexType index;
    FloatImageType::IndexType neighbourindex;

    HessianFilter::OutputImageType::PixelType current_hessian;


    std::cout << "Starting calculating non Maximum Supression" << std::endl;

    x_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(0);
    y_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(1);
    z_size = maxMedialnessOverScales->GetLargestPossibleRegion().GetSize(2);

    float a[3][3];
    double* a_double;
    float b[3][3];
    float d[3];

    float *float_ptr;
    double *double_ptr;

    float current_value;

    for(x = 1; x < x_size - 1; x++)
    {
      index[0] = x;

      if(x % 10 == 0)
        std::cout << "x " << x << std::endl;

      for(y = 1; y < y_size - 1; y++)
      {
        index[1] = y;

        for(z = 1; z < z_size - 1; z++)
        {
          index[2] = z;

          current_value = maxMedialnessOverScales->GetPixel(index);

          current_hessian = hessianImages.back()->GetPixel(index);


          //current_hessian.FixedArray(a_double);
          a_double = current_hessian.GetDataPointer();

          int u=0;
          for(float_ptr = &a[0][0], double_ptr = a_double, u=0; u<6; u++, float_ptr++, double_ptr++)
              *float_ptr = (float)*double_ptr;

          NumericsHelper::calculateEigenDecompositionSymmetric3x3(a,b,d);



        }
      }
    }
}

*/



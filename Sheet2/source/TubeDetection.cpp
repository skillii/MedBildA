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



#define START_LEVEL (1) //number of scale-levels to skip

const int TubeDetection::scale_levels = SCALE_LEVELS;
const float TubeDetection::tube_r[SCALE_LEVELS] = {1.5, 1.5, 1.5};
const int TubeDetection::alpha_steps = 16;



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

        /*
        gradX = imagePyramid.at(i)->Clone();
        gradY = imagePyramid.at(i)->Clone();
        gradZ = imagePyramid.at(i)->Clone();

        gradX->Allocate();
        gradY->Allocate();
        gradZ->Allocate();*/

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

                    /*
                    gradX->SetPixel(index_f, gradient.GetDataPointer()[0]);
                    gradY->SetPixel(index_f, gradient.GetDataPointer()[1]);
                    gradZ->SetPixel(index_f, gradient.GetDataPointer()[2]); */

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

void TubeDetection::calcMedialness()
{
    unsigned x,y,z;
    HessianFilter::OutputImageType::IndexType index;
    unsigned int x_size;
    unsigned int y_size;
    unsigned int z_size;

    HessianFilter::OutputImagePixelType current_hessian;

    float a[3][3];
    double* a_double;
    float b[3][3];
    float d[3];

    float *float_ptr;
    double *double_ptr;
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

                    current_hessian = hessianImages.at(i)->GetPixel(index);


                    //current_hessian.FixedArray(a_double);
                    a_double = current_hessian.GetDataPointer();

                    int u=0;
                    for(float_ptr = &a[0][0], double_ptr = a_double, u=0; u<6; u++, float_ptr++, double_ptr++)
                        *float_ptr = (float)*double_ptr;

                    NumericsHelper::calculateEigenDecompositionSymmetric3x3(a,b,d);

                    float medialness = calcMedialness(i, x, y, z, d, b);

                    medialnessimage->SetPixel(index, medialness);
                }
            }
        }

        medialnessImages.push_back(medialnessimage);
    }

    std::cout << std::endl << "Finished calculation of medialness" << std::endl;
}

float TubeDetection::calcMedialness(unsigned level, unsigned x, unsigned y, unsigned z, float ew[3], float ev[3][3])
{
    int i,j;


    float alpha_step = 2*M_PI/alpha_steps;
    float alpha = 0;

    float c[alpha_steps];

    FloatImageType::Pointer img = imagePyramid.at(level);

    itk::Vector<float, 3> v, grad;
    itk::Vector<float, 3> v_i[2];



    int min_ev;
    NumericsHelper::min(ew, 3, &min_ev);


    for(i = 0,j = 0; i < 3; i++)
    {
        if(i != min_ev)
        {
            v_i[j].SetElement(0, ev[i][0]);
            v_i[j].SetElement(1, ev[i][1]);
            v_i[j].SetElement(2, ev[i][2]);
            j++;
        }
    }




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



    float final = s-rc;

    if(final < 0)
        final = 0;

    return final;
}





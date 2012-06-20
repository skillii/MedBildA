/*
 * TubeDetection.h
 *
 *  Created on: Jun 16, 2012
 *      Author: tom
 */

#ifndef TUBEDETECTION_H_
#define TUBEDETECTION_H_

#include "Definitions.h"
#include <vector>

#define SCALE_LEVELS 3


class TubeDetection {
    FloatImageType::Pointer lungSegment;
    FloatImageType::Pointer raw;
    FloatImageType::Pointer lung;
    FloatImageType::Pointer maxMedialnessOverScales;


    //Eigenvector with Smallest eigenvalue (x = eigenvector[0],...)
    FloatImageType::Pointer eigenvector[3];


    std::vector<FloatImageType::Pointer> imagePyramid;
    std::vector<FloatImageType::Pointer> medialnessImages;
    std::vector<HessianFilter::OutputImageType::Pointer> hessianImages;
    //std::vector<GradientImageFilter::OutputImageType::Pointer> gradientImages;

    std::vector<FloatImageType::Pointer> gradientX, gradientY, gradientZ;




    static const int scale_levels;
    static const float tube_r[SCALE_LEVELS];



    void allocateEigenvectorImage();
    float calcMedialness(unsigned level, unsigned x, unsigned y, unsigned z, float ew[3], float ev[3][3]);
    void eigenValueDecomposition(int scale, HessianFilter::OutputImageType::IndexType index, float V[3][3], float d[3]);

public:

    static const int alpha_steps;

    void cropLung();
    void buildImagePyramid();
    void calcHessian();
    void calcGradients();

    void calcMedialness();
    void calcMaxMedialness();

    TubeDetection(FloatImageType::Pointer lungSegment, FloatImageType::Pointer img);
    virtual ~TubeDetection();

    std::vector<FloatImageType::Pointer> getImgPyramid()
    {
        return imagePyramid;
    }

    std::vector<FloatImageType::Pointer> getMedialnessImages()
    {
        return medialnessImages;
    }

    FloatImageType::Pointer getMaxMedialnessImage()
    {
        return this->maxMedialnessOverScales;
    }

    std::vector<FloatImageType::Pointer> getGradientImage(int dim)
    {
      switch(dim)
      {
        case 0: return gradientX;
        case 1: return gradientY;
        case 2: return gradientZ;
        default: throw 0;
      }
    }


    FloatImageType::Pointer getEVImage(int dim)
    {
      return eigenvector[dim];
    }


};

#endif /* TUBEDETECTION_H_ */

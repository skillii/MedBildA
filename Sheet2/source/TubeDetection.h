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

#define SCALE_LEVELS 5


class TubeDetection {
    FloatImageType::Pointer lungSegment;
    FloatImageType::Pointer raw;
    FloatImageType::Pointer lung;

    std::vector<FloatImageType::Pointer> imagePyramid;
    std::vector<FloatImageType::Pointer> medialnessImages;
    std::vector<HessianFilter::OutputImageType::Pointer> hessianImages;
    //std::vector<GradientImageFilter::OutputImageType::Pointer> gradientImages;

    std::vector<FloatImageType::Pointer> gradientX, gradientY, gradientZ;


    static const int scale_levels;
    static const float tube_r[SCALE_LEVELS];
    static const int alpha_steps;



    float calcMedialness(unsigned level, unsigned x, unsigned y, unsigned z, float ew[3], float ev[3][3]);

public:

    void cropLung();
    void buildImagePyramid();
    void calcHessian();
    void calcGradients();

    void calcMedialness();

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

    std::vector<FloatImageType::Pointer> getGradientImage(int dim)
    {
      switch(dim)
      {
        case 1: return gradientX;
        case 2: return gradientY;
        case 3: return gradientZ;
        default: throw 0;
      }
    }
};

#endif /* TUBEDETECTION_H_ */

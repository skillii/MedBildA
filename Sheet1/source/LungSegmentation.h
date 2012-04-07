/*
 * LungSegmentation.h
 *
 *  Created on: Apr 4, 2012
 *      Author: tom
 */

#ifndef LUNGSEGMENTATION_H_
#define LUNGSEGMENTATION_H_

#include "Definitions.h"




/**
 * This class applies a thresholding filter, uses shape labeling to find
 * connected components and finally chooses lung according to shape sizes.
 */
class LungSegmentation {
	FloatImageType::Pointer img_;
	FloatImageType::Pointer thresholdimg_;
	FloatImageType::Pointer lungimg_;


	/**
	 * Removes all shape labels, except the lung
	 *
	 * @param shape_labels the shape labels
	 */
	void SelectLungShape(ConnectedComponentFilter::OutputImageType* shape_labels);


public:
	LungSegmentation(FloatImageType::Pointer img);
	virtual ~LungSegmentation();


	/**
	 * performs the lung segmentation
	 */
	void run();

	/**
	 * returns the threshold image(binary image), which is an intermediate result
	 * of the lung segmentation
	 *
	 * @return the threshold image (binary image)
	 */
	FloatImageType::Pointer getThresholdImage()
	{
		return thresholdimg_;
	}


	/**
	 * returns the image containing the lung
	 *
	 * @return image containing the lung
	 */
	FloatImageType::Pointer getLungImage()
	{
		return lungimg_;
	}
};

#endif /* LUNGSEGMENTATION_H_ */

/*
 * LungSegmentation.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: tom
 */

#include <cmath>
#include <iostream>


#include "Definitions.h"
#include "LungSegmentation.h"




LungSegmentation::LungSegmentation(FloatImageType::Pointer img)
{
	this->img_ = img;
	this->thresholdimg_ = NULL;
	this->lungimg_ = NULL;
}

LungSegmentation::~LungSegmentation()
{
}


void LungSegmentation::run()
{
	//==================================================================
	//Perform Thresholding
	//==================================================================
	ThresholdFilterType::Pointer t_filter = ThresholdFilterType::New();

    t_filter->SetInput(this->img_);

	t_filter->SetOutsideValue(255);
	t_filter->SetInsideValue(0);

	t_filter->SetNumberOfHistogramBins(128);

	t_filter->Update();

	this->thresholdimg_ = t_filter->GetOutput();



	//==================================================================
	// Find the connected components of the BinaryImage(Thresholdimage):
	//==================================================================

	ConnectedComponentFilter::Pointer c_filter = ConnectedComponentFilter::New();

	c_filter->SetInput(this->thresholdimg_);
	c_filter->SetInputForegroundValue(0);
	c_filter->Update();

	ConnectedComponentFilter::OutputImageType::Pointer shape_labels = c_filter->GetOutput();

	//==================================================================
	// Extract the lung (using shape size as classifier)
	//=================================================================

	SelectLungShape(shape_labels);

	//==================================================================
	// Convert the label map to an image
	//=================================================================

	LabelMapToLabelImageFilter::Pointer labelmaptoimage_filter = LabelMapToLabelImageFilter::New();

	labelmaptoimage_filter->SetInput(shape_labels);
	labelmaptoimage_filter->Update();

	this->lungimg_ = labelmaptoimage_filter->GetOutput();
}


void LungSegmentation::SelectLungShape(ConnectedComponentFilter::OutputImageType* shape_labels)
{
	std::cout << "number of labels:" << shape_labels->GetNumberOfLabelObjects() << std::endl;

	if(shape_labels->GetNumberOfLabelObjects() < 2)
		return;

	//look for the second largest label, since this is the lung!

	ConnectedComponentFilter::OutputImageType::LabelObjectType::Pointer max_label = NULL;
	ConnectedComponentFilter::OutputImageType::LabelObjectType::Pointer max2_label = NULL;
	unsigned max = 0;
	unsigned max2 = 0;

	for(unsigned i = 0; i < shape_labels->GetNumberOfLabelObjects(); i++)
	{
		ConnectedComponentFilter::OutputImageType::LabelObjectType::Pointer label = shape_labels->GetNthLabelObject(i);

		unsigned size = label->GetNumberOfPixels();

		if(size >= max)
		{
			max2_label = max_label;
			max2 = max;

			max_label = label;
			max = size;
		}
		else if(size >= max2)
		{
			max2_label = label;
			max2 = size;
		}

		std::cout << "size of " << i << "th label:" << size << std::endl;
	}

	//now remove all labels, exept the second largest one:

	unsigned offset = 0;

	while(shape_labels->GetNumberOfLabelObjects() > 1)
	{
		ConnectedComponentFilter::OutputImageType::LabelObjectType::Pointer label = shape_labels->GetNthLabelObject(offset);

		if(label->GetLabel() == max2_label->GetLabel())
		{
			offset++;
		}
		else
		{
			shape_labels->RemoveLabel(label->GetLabel());
		}
	}
}



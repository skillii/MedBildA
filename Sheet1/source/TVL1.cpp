/*
 * TVL1.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: tom
 */


#include "Definitions.h"
#include "TVL1.h"
#include <cmath>



int TVL1::nrIterations = 500;
float TVL1::lambda = 0.3;
float TVL1::tau = 0.02;
float TVL1::sigma;

TVL1::TVL1(FloatImageType::Pointer img) {
	 TVL1::sigma = 1/sqrt(12*tau);
	 this->img = img;
}

TVL1::~TVL1() {
}


void TVL1::Denoise(void)
{
	unsigned x,y,z;

	FloatImageType::SizeType size;
	PrimalImageType::SizeType size_p;

	PrimalImageType::Pointer p = PrimalImageType::New();
	FloatImageType::Pointer u;
	FloatImageType::Pointer u_quer = FloatImageType::New();

	size = img->GetLargestPossibleRegion().GetSize();
	size_p[0] = size[0]; size_p[1] = size[1]; size_p[2] = size[2];
	size_p[3] = 3;	//vector of length 3 for each pixel in the image


	p->SetRegions(size_p);
	p->Allocate();


	//u->SetRegions(size);
	//u->Allocate();
	u = this->img->Clone();

	u_quer->SetRegions(size);
	u_quer->Allocate();

	int iteration;

	for(iteration = 0; iteration < nrIterations; iteration++)
	{
		for(x = 0; x < img->GetLargestPossibleRegion().GetSize(0); x++)
		{
			for(y = 0; y < img->GetLargestPossibleRegion().GetSize(1); y++)
			{
				for(z = 0; z < img->GetLargestPossibleRegion().GetSize(2); z++)
				{
					FloatImageType::IndexType index;
					index[0] = x; index[1] = y; index[2] = z;


				}
			}
		}
	}
}


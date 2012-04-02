/*
 * TVL1.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: tom
 */


#include "Definitions.h"
#include "TVL1.h"
#include <cmath>



static int TVL1::nrIterations = 500;
static float TVL1::lambda = 0.3;
static float TVL1::tau = 0.02;
static float TVL1::sigma;

TVL1::TVL1(FloatImageType::Pointer img) {
	 TVL1::sigma = 1/sqrt(12*tau);
	 this->img = img;
}

TVL1::~TVL1() {
}


void TVL1::Denoise(void)
{
	int x,y,z;

	FloatImageType::IndexType start;
	FloatImageType::SizeType size;

	FloatImageType::Pointer p = FloatImageType::New();
	FloatImageType::Pointer u = FloatImageType::New();
	FloatImageType::Pointer u_quer = FloatImageType::New();

	start[0] = 0; start[1] = 0; start[2] = 0;
	size = img->GetLargestPossibleRegion().GetSize();

	FloatImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);


	p->SetRegion(region);
	p->Allocate();


	u->SetRegion(region);
	u->Allocate();

	u_quer->SetRegion(region);
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
					FloatImageType:IndexType index;
					index[0] = x; index[1] = y; index[2] = z;


				}
			}
		}
	}
}


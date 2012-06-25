/*
 * TVL1.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: tom
 */


#include "Definitions.h"
#include "TVL1.h"
#include <cmath>
#include <iostream>




int TVL1::nrIterations = 500;
float TVL1::lambda = 0.3;
float TVL1::tau = 0.02;
float TVL1::sigma;
float TVL1::tau_times_lambda;
float TVL1::theta = 0.7; //Overrelaxation constant

TVL1::TVL1(FloatImageType::Pointer img, float max_value)
{
	 TVL1::sigma = 1/sqrt(12*tau);
	 TVL1::tau_times_lambda = TVL1::tau * TVL1::lambda;
	 this->img = img;
	 this->max_value = max_value;
}

TVL1::~TVL1()

{
}


FloatImageType::Pointer TVL1::Denoise(void)
{
	unsigned x,y,z;

	FloatImageType::SizeType size;
	PrimalImageType::SizeType size_p;
	PrimalImageType::IndexType index_p;

	PrimalImageType::Pointer p = PrimalImageType::New();

	size = img->GetLargestPossibleRegion().GetSize();
	size_p[0] = size[0]; size_p[1] = size[1]; size_p[2] = size[2];
	size_p[3] = 3;	//vector of length 3 for each pixel in the image


	index_p[0] = 0; index_p[1] = 0; index_p[2] = 0; index_p[3] = 0;


	//==================================================================
	//BUILD PRIMAL IMAGE, Start value = 0
	//==================================================================
	p->SetRegions(size_p);
	p->Allocate();

    //Set pixels in p image to 0
	//PrimalImageType::PixelType pixelValue = 0.0f;
	p->FillBuffer(0);


	//===================================================================
	//BUILD DUAL IMAGES, Start values = f
	//===================================================================
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(this->img);
    duplicator->Update();

    FloatImageType::Pointer u      = duplicator->GetOutput();

    duplicator = DuplicatorType::New();
	duplicator->SetInputImage(this->img);
    duplicator->Update();
    FloatImageType::Pointer u_dash = duplicator->GetOutput();


    FloatImageType::IndexType index;
    unsigned int x_size = img->GetLargestPossibleRegion().GetSize(0);
    unsigned int y_size = img->GetLargestPossibleRegion().GetSize(1);
    unsigned int z_size = img->GetLargestPossibleRegion().GetSize(2);



    //TODO ?????????????
    //u->Allocate();
   // u_dash->Allocate();


   //===================================================================
    //ITERATIVE PRIMAL DUAL UPDATES
    //===================================================================

    int iteration;
    FloatImageType::PixelType u_value;
    FloatImageType::PixelType f_value;
    FloatImageType::PixelType dx;
    FloatImageType::PixelType dy;
    FloatImageType::PixelType dz;


    FloatImageType::IndexType index_dx;
    FloatImageType::IndexType index_dy;
    FloatImageType::IndexType index_dz;


    PrimalImageType::IndexType index_dx_bw;
    PrimalImageType::IndexType index_dy_bw;
    PrimalImageType::IndexType index_dz_bw;

    float p_temp[3];
    float p_norm;
    float u_temp;

    index_dx_bw[3] = 0;
    index_dy_bw[3] = 1;
    index_dz_bw[3] = 2;

	for(iteration = 0; iteration < nrIterations; iteration++)
	{
		std::cout << "Iteration #" << iteration << std::endl;

		for(x = 0; x < x_size ; x++)
		{
			index[0]    = x;
			index_dx[0] = x+1;
			index_dy[0] = x;
			index_dz[0] = x;
			index_p[0] = x;

			index_dx_bw[0] = x - 1;
			index_dy_bw[0] = x;
			index_dz_bw[0] = x;


			for(y = 0; y < y_size; y++)
			{
				index[1] = y;
				index_dx[1] = y;
				index_dy[1] = y+1;
				index_dz[1] = y;
				index_p[1] = y;

				index_dx_bw[1] = y;
				index_dy_bw[1] = y - 1;
				index_dz_bw[1] = y;

			    for(z = 0; z < z_size; z++)
				{
			    	 //std::cout << "Coordinate " << x << "/" << y << "/" << z << std::endl;
				    //Indices for gradient, forward differences
				    index[2] = z;
				    index_dx[2] = z;
				    index_dy[2] = z;
				    index_dz[2] = z+1;
				    index_p[2] = z;

				    //=========================================================
				    //Calculate gradient of u_dash
				    //=========================================================
				    dx = dy = dz = 0;

				    if(x+1 < x_size)
				      dx = u_dash->GetPixel(index_dx) -
				           u_dash->GetPixel(index);


				    if(y+1 < y_size)
				      dy = u_dash->GetPixel(index_dy) -
				           u_dash->GetPixel(index);

				    if(z+1 < z_size)
				      dz = u_dash->GetPixel(index_dz) -
				           u_dash->GetPixel(index);

                    dx *= sigma; dy *= sigma; dz *= sigma;

                    //std::cout << "----Gradient calculated" << std::endl;

                    //=========================================================
                    //Calculate dual update
                    //=========================================================
                    index_p[3] = 0;
                    p_temp[0] = p->GetPixel(index_p) + dx;

                    index_p[3] = 1;
                    p_temp[1] = p->GetPixel(index_p) + dy;

                    index_p[3] = 2;
                    p_temp[2] = p->GetPixel(index_p) + dz;

                    //std::cout << "----Dual update calculated" << std::endl;

                    //=========================================================
                    //Project p and update p image
                    //=========================================================
                    p_norm = sqrt(p_temp[0]*p_temp[0] + p_temp[1]*p_temp[1] + p_temp[2]*p_temp[2]);

                    if(p_norm < 1)
                    {
                      index_p[3] = 0;  p->SetPixel(index_p, p_temp[0]);
                      index_p[3] = 1;  p->SetPixel(index_p, p_temp[1]);
                      index_p[3] = 2;  p->SetPixel(index_p, p_temp[2]);
                    }
                    else
                    {
                      index_p[3] = 0;  p->SetPixel(index_p, p_temp[0] / p_norm);
                      index_p[3] = 1;  p->SetPixel(index_p, p_temp[1] / p_norm);
                      index_p[3] = 2;  p->SetPixel(index_p, p_temp[2] / p_norm);
                    }

                    //std::cout << "----P projected and updated" << std::endl;
                    //=========================================================
                    //Calculate primal update (u_temp)
                    //=========================================================
                    index_p[3] = 0;
                    index_dx_bw[2] = z;
                    index_dy_bw[2] = z;
                    index_dz_bw[2] = z - 1;


                    if(x - 1 > 0)
                      u_temp = p->GetPixel(index_p) - p->GetPixel(index_dx_bw);
                    else
                      u_temp = p->GetPixel(index_p);

                    index_p[3] = 1;

                    if(y - 1 > 0)
                      u_temp += (p->GetPixel(index_p) - p->GetPixel(index_dy_bw));
                    else
                      u_temp += p->GetPixel(index_p);

                    index_p[3] = 2;

                    if(z - 1 > 0)
                      u_temp += (p->GetPixel(index_p) - p->GetPixel(index_dz_bw));
                    else
                      u_temp += p->GetPixel(index_p);

                    u_temp *= tau;
                    u_temp += u->GetPixel(index);

                    //std::cout << "----Primal update calculated" << std::endl;
                    //=========================================================
                    //Perform resolvent step for u
                    //=========================================================
                    u_value = u->GetPixel(index);
                    f_value = this->img->GetPixel(index);


                    #ifdef USE_PRIMAL_DUAL_ROF
                      u->SetPixel(index, (u_temp + tau*lambda*f_value) / (1+ tau*lambda));

                    #endif
                    /* TVL1 Variante 1...

                    if(u_value - f_value > TVL1::tau_times_lambda)
                      u->SetPixel(index, u_value - TVL1::tau_times_lambda);

                    if(u_value - f_value < -TVL1::tau_times_lambda)
                      u->SetPixel(index, u_value + TVL1::tau_times_lambda);

                    if(abs(u_value - f_value) <= TVL1::tau_times_lambda)
                      u->SetPixel(index, f_value); */


                    #ifdef USE_TVL1

                    if(u_temp - f_value > TVL1::tau_times_lambda)
                      u->SetPixel(index, u_temp - TVL1::tau_times_lambda);

                    else if(u_temp - f_value < -TVL1::tau_times_lambda)
                      u->SetPixel(index, u_temp + TVL1::tau_times_lambda);

                    else //if(abs(u_temp - f_value) <= TVL1::tau_times_lambda)
                      u->SetPixel(index, f_value);

                    #endif

                    //std::cout << "----Resolvent step calculated" << std::endl;

                    //=========================================================
                    //Finally perform overrelaxation
                    //=========================================================
                    u_dash->SetPixel(index,u->GetPixel(index) + theta*(u->GetPixel(index) - u_value));

                    //std::cout << "----Overrelaxation calculated" << std::endl;
				}
			}
		}
	}

	duplicator = DuplicatorType::New();
	duplicator->SetInputImage(u_dash);
    duplicator->Update();

    std::cout << "Done with calculations!" << std::endl;

    this->img = duplicator->GetOutput();

    return this->img;
}


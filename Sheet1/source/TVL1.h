/*
 * TVL1.h
 *
 *  Created on: Apr 2, 2012
 *      Author: tom
 */

#ifndef TVL1_H_
#define TVL1_H_

class TVL1 {
	FloatImageType::Pointer img;
    float max_value;


	static int nrIterations;
	static float lambda, tau, sigma;
	static float tau_times_lambda;
	static float theta;
public:
	TVL1(FloatImageType::Pointer img, float max_value);
	virtual ~TVL1();

	FloatImageType::Pointer Denoise();
};

#endif /* TVL1_H_ */

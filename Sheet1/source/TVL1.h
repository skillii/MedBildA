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



	static int nrIterations;
	static float lambda, tau, sigma;
public:
	TVL1(FloatImageType::Pointer img);
	virtual ~TVL1();

	void Denoise();
};

#endif /* TVL1_H_ */

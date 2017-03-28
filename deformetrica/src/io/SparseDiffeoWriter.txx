/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoWriter_txx
#define _SparseDiffeoWriter_txx

#include "SparseDiffeoWriter.h"

#include "itksys/SystemTools.hxx"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
SparseDiffeoWriter<TScalar, Dimension>
::SparseDiffeoWriter()
 {

 }



template <class TScalar, unsigned int Dimension>
SparseDiffeoWriter<TScalar, Dimension>
::~SparseDiffeoWriter()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
SparseDiffeoWriter<TScalar, Dimension>
::SetKernelType(const char* type)
 {
	if (itksys::SystemTools::Strucmp(type,"Exact") == 0)
		m_KernelType = Exact;
	else if (itksys::SystemTools::Strucmp(type,"P3M") == 0)
		m_KernelType = P3M;
	else if (itksys::SystemTools::Strucmp(type,"FGT") == 0)
		m_KernelType = FGT;
	else
		throw std::runtime_error("unknown kernel type");

	return;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
SparseDiffeoWriter<TScalar, Dimension>
::Update()
 {

	std::ofstream outfile(m_FileName);

	outfile << "Sparse Diffeo file v1.0" << std::endl;

	outfile << "Kernel Width = " << std::endl << m_Def->GetKernelWidth() << std::endl;
	// outfile << "Kernel Type = " << std::endl << m_KernelType << std::endl;
	outfile << "Number of Time Points = " << std::endl << m_Def->GetNumberOfTimePoints() << std::endl;
	outfile << "DataDomain = " << std::endl << m_Def->GetDataDomain() << std::endl;
	outfile << "Padding Factor = " << std::endl << m_Def->GetPaddingFactor() << std::endl;

	m_Def->Update();
	if (m_Def->OutOfBox())
		std::cerr << "Warning: you are computing an corrupted diffeos (outofbox). Results may be inaccurate!" << std::endl;

	MatrixList XT = m_Def->GetTrajectoryPositions();
	MatrixList MomT = m_Def->GetTrajectoryMomentas();

	int numCP = XT[0].rows();
	outfile << "Number of Points = " << std::endl << numCP << std::endl;

	outfile << "Positions = " << std::endl;
	for (int t = 0; t < m_Def->GetNumberOfTimePoints(); t++)
	{
		outfile << "t = " << t << std::endl;
		for (int i = 0; i < numCP; i++)
			outfile << XT[t].get_row(i) << std::endl;
	}


	outfile << "Momentas = " << std::endl;
	for (int t = 0; t < m_Def->GetNumberOfTimePoints(); t++)
	{
		outfile << "t = " << t << std::endl;
		for (int i = 0; i < numCP; i++)
			outfile << MomT[t].get_row(i) << std::endl;
	}

	outfile.close();

	return;

 }



#endif /* _SparseDiffeoWriter_txx */

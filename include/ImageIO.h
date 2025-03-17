#pragma once
#include <string>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <bitset>

using namespace std;

template <typename T> static T** AllocateImage(int ncols, int nrows);
template <typename T> static void FreeImage(T** image);

template <typename T> T** AllocateImage(int ncols, int nrows)
{
	try
	{
		T** ptr = new T * [nrows];
		T* pool = new T[nrows * ncols];
		for (unsigned i = 0; i < nrows; ++i, pool += ncols)
			ptr[i] = pool;
		return ptr;
	}
	catch (std::bad_alloc)
	{
		std::cerr << "Not enough memory could be allocated by the system." << std::endl;    
		assert(false);
		return NULL;
	}
}

template <typename T> void FreeImage(T** image)
{
	if (image == NULL)
		return;

	if (image != NULL)
	{
		delete[] image[0];
		delete[] image;
		image = NULL;
	}
	delete[] image;
	image = NULL;
}


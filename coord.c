#include <assert.h>
#include <stdlib.h>
#include "coord.h"

size_t golden_spiral_distribution(coord * vectors, size_t num_vectors)
{
	// sanity checks
	assert(vectors != NULL);

	// must be an odd number
	if(num_vectors % 2 == 0)
		num_vectors--;

	// calculate the sampling vectors
	int kmax = (int) (num_vectors - 1) / 2;
	double part1, part2, part3;
	for(int k = -kmax, i=0; k <= kmax; k++, i++)
	{
		part1 = 2.0 * k / num_vectors;
		part2 = 2.0 * k * coord_pi / coord_golden_ratio;
		part3 = cos(asin(part1));
		vectors[i].x = part3 * cos(part2);
		vectors[i].y = part3 * sin(part2);
		vectors[i].z = part1;
	}

	return num_vectors;
}
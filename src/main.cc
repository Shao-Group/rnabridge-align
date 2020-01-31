/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "config.h"
#include "previewer.h"
#include "coral.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));

	if(argc == 1)
	{
		print_copyright();
		print_help();
		printf("\n");
		return 0;
	}

	parse_arguments(argc, argv);

	if(verbose >= 1)
	{
		print_copyright();
		printf("\n");
		print_command_line(argc, argv);
		printf("\n");
	}

	previewer pv;
	pv.preview();

	if(preview_only == true) return 0;

	coral cr;
	cr.resolve();

	return 0;
}

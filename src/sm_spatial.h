/* Copyright (C) 2009-2012 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/


#ifndef SM_SPATIAL_H_
#define SM_SPATIAL_H_




class Stringmol_Spatial: public stringPM
{

private:

public:

	int Run(int argc, char *argv[]);

	int ConfigureFromFile(const char *fn, stringPM *A, smsprun **run, int runno);
	int GridSavePNG(stringPM *A, smpic pt);

};


#endif //SM_SPATIAL_H_

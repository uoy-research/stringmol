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


#ifndef OPCODES_H_
#define OPCODES_H_


int OpcodeCopy(s_ag *act, const bool domut,float indelrate,
		float subrate, const unsigned int maxl,
		swt	*blosum, const int granular_1, long &biomass
		);

int OpcodeComassCopy(s_ag *act, const bool domut,float indelrate,
		float subrate, const unsigned int maxl,
		swt	*blosum, const int granular_1, long &biomass,
		int *mass);

#endif

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


char * OpcodeSearchInner(char *iptr, char *sp, swt *T, const int *itog,
		int *ftog, int maxl);

void OpcodeSearch(s_ag *act, swt *blosum, const unsigned short int maxl);

void OpcodeMove(s_ag *act);

void OpcodeIncrementRead(s_ag *act, bool granular_1);

void OpcodeToggle(s_ag *act);

int OpcodeCleave(s_ag *act, s_ag *nexthead, SMspp *spl,
		unsigned long int *agct,
		const unsigned int timestep, const unsigned int maxl0);

char * OpcodeIf(char *ip, char *rp, char *sp, swt *T, const int maxl);

void OpcodeInsertInstruction(const s_ag * act, int inst_idx,
		int *mass, swt * blosum,
		const int writePtrOpcodeIndex = -1);

int OpcodeCopy(      s_ag *act, const bool domut,float indelrate,
		float subrate, const unsigned int maxl,
		swt	*blosum, const int granular_1, long &biomass
		);

int OpcodeComassCopy(s_ag *act, const bool domut,float indelrate,
		float subrate, const unsigned int maxl,
		swt	*blosum, const int granular_1, long &biomass,
		int *mass);



#endif /* OPCODES_H_ */

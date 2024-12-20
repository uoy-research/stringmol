/* Copyright (C) 2009-2012 Verena Fischer and Simon Hickinbotham        */
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


#ifdef __cplusplus
extern "C" {
#endif


#ifndef SNIPPET_H_
#define SNIPPET_H_


int ga_step_int(int **pop, const double *eval, const int POPSIZE, const int PARAMETERS, const int minval, const int maxval, int *wn);


#endif /* SNIPPET_H_ */

#ifdef __cplusplus
}
#endif

/* 
 * TwinBlad
 *
 * Copyright (c) 2014 - 2015 Patrick Rooney (darraghrooney@gmail.com)
 *
 * This file is part of TwinBlad.
 *
 * TwinBlad is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TwinBlad is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TwinBlad.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef MENUS_INCLUDED
#define MENUS_INCLUDED

#include "LSystem2d.h"

// mainMenu is run automatically
char mainMenu(void);

// From main menu, user chooses one of five secondary menus
void randomMenu(void);
void PauliMenu(void);
void realdiagMenu(void);
void pmzMenu(void);
void fileMenu(void);

#endif // MENUS_INCLUDED

/***************************************************************************
 *   Copyright (C) 2009 by Hermann Sonntag   *
 *   hermann.sonntag@tu-ilmenau.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
struct argos_args_info {
	size_t STORED_SAMPLING_RATE,
	       size_t	STORED_PKT_SIZE,
	       size_t MEAS_NO_OF_PACKETS,
	       size_t STORED_MAG_CHANNELS,
	       size_t STORED_MAG_REF_CHANNELS,
	       size_t STORED_ELEC_CHANNELS,
	       size_t MAG_SENSOR_A_POSITION,
	       size_t MAG_SENSOR_A_ORIENTATION,
	       size_t MAG_SENSOR_B_POSITION,
	       size_t MAG_SENSOR_B_ORIENTATION,
	       size_t MAG_SENSOR_C_POSITION,
	       size_t MAG_SENSOR_C_ORIENTATION,
	       size_t MEAS_START_PACKET_NO 
}
int argos_read_header(struct argos_args_info * args_info, const char * filename);

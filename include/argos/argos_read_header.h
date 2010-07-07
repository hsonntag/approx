/*  argos/argos_read_header.h                                              *
 ***************************************************************************
 *   Copyright (C) 2009, 2010 Hermann Sonntag                              *
 *   hermann.sonntag@tu-ilmenau.de                                         *
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
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 ***************************************************************************/

struct argos_args_info {
	size_t stored_sampling_rate,
	        stored_pkt_size,
	        meas_no_of_packets,
	        stored_mag_channels,
	        stored_mag_ref_channels,
	        stored_elec_channels,
	        mag_sensor_a_position,
	        mag_sensor_a_orientation,
	        mag_sensor_b_position,
	        mag_sensor_b_orientation,
	        mag_sensor_c_position,
	        mag_sensor_c_orientation,
	        meas_start_packet_no 
}

int argos_read_header(struct argos_args_info * args_info, const char * filename);


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
	size_t stored_sampling_rate,
	       size_t	stored_pkt_size,
	       size_t meas_no_of_packets,
	       size_t stored_mag_channels,
	       size_t stored_mag_ref_channels,
	       size_t stored_elec_channels,
	       size_t mag_sensor_a_position,
	       size_t mag_sensor_a_orientation,
	       size_t mag_sensor_b_position,
	       size_t mag_sensor_b_orientation,
	       size_t mag_sensor_c_position,
	       size_t mag_sensor_c_orientation,
	       size_t meas_start_packet_no 
}
int argos_read_header(struct argos_args_info * args_info, const char * filename);

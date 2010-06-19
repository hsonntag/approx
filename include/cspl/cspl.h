/*   cspl/cspl.h                                                           *
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
#ifndef __CSPL_H__
#define __CSPL_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <cspl/cspl_qrs.h>
#include <cspl/cspl_qrs_fit.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS /* Hier beginnen die Deklarationen für cspl. */

/**
 * Diese Methode berechnet die Fouriertransformation eines reellen Datensatzes
 * und speichert diese im halbkomplexen Format.
 * @author Hermann Sonntag
 * @param data Pointer zu dem Datensatz der transformiert wird,
 * wobei die Transformation an gleicher Stelle abgelegt wird
 * @param size Länge des zu transformierenden Datensatzes
 * @return Status der Transformation als int
 */
int cspl_real_fft (double * data, size_t size);

/**
 * Diese Methode berechnet die Absolutwerte aus einem Datensatz im halbkomplexen Format.
 * @author Hermann Sonntag
 * @param data Pointer zu einem halbkomplexen Datensatz
 * @param size Länge des halbkomplexen Datensatzes
 * @return Status der Berechnung
 */
int cspl_halfcomplex_abs (double * data, size_t size);

/**
 * Diese Methode berechnet die Summe der Differnzenquadrate zweier diskreter Funktionen ungleicher Länge
 * über die Länge der zweiten Funktion, zu jeder Verschiebung der Funktionen zueinander.
 * @author Hermann Sonntag
 * @param rms Ergebnisvektor
 * @param x erste diskrete Funktion
 * @param y zweite diskrete Funktion
 * @param length Länge jeder einzelnen Summenberechnung, welche der Länge der zweiten Funktion entspricht
 * @param size Länge der ersten Funktion
 * @return Status der Berechnung
 */
int cspl_root_mean_square (double * rms, double * x, double * y, size_t length, size_t size);

/**
 * Diese Methode normiert eine diskrete Funktion und macht diese mittelwertfrei.
 * @param signal Pointer zu der Funktion, welche an dieser Stelle normiert wird.
 * @param size Länge der Funktion
 * @return Skalierungsdivident der Normierung
 */
int cspl_norm (double * signal, size_t size);

/**
 * Diese Methode erstellt die Liste aller lokalen Minima in einem bestimmten Bereich,
 * welche unter einem Referenzminimum liegen.
 * Wenn zwei lokale Minima nur eine Bereichslänge voneinander entfernt gefunden werden,
 * wird die Suche abgebrochen.
 * @author Hermann Sonntag
 * @param min_n Liste, in der die Minima gespeichert werden
 * @param signal diskrete Funktion, in der die Minima gesucht werden
 * @param size Länge der Funktion
 * @param length Bereichslänge
 * @param min Referenzminimum
 * @return Anzahl der gefundenen Mimima
 */
int cspl_eval_periodic_min2 (unsigned int * min_n, double * signal, size_t size, size_t length, double min); 

/**
 * Diese Methode erstellt die Liste aller lokalen Maxima in einem bestimmten Bereich,
 * welche unter einem Referenzmaximum liegen.
 * Wenn zwei lokale Maxima nur eine Bereichslänge voneinander entfernt gefunden werden,
 * wird die Suche abgebrochen.
 * @author Hermann Sonntag
 * @param min_n Liste, in der die Maxima gespeichert werden
 * @param signal diskrete Funktion, in der die Maxima gesucht werden
 * @param size Länge der Funktion
 * @param length Bereichslänge
 * @param min Referenzmaximum
 * @return Anzahl der gefundenen Maxima
 */
int cspl_eval_periodic_max2 (unsigned int * max_n, double * signal, size_t size, size_t length, double max); 

int cspl_radix2_xcorr(double * c_xy, double * x, double * y, size_t size);

/**
 * Diese Methode berechnet die Kreuzkorrelation zweier Funkionen unterschiedlicher Länge.
 * @author Hermann Sonntag
 * @param c_xy Pointer zu dem Datensatz, in dem die Kreuzkorrelation gespeichert wird
 * @param x erste Funktion
 * @param y zweite Funktion
 * @param length Länge der zweiten Funktion
 * @param size Länge der ersten Funktion
 * @return Status der Korrelationsberechnung
 */
int cspl_xcorr(double * c_xy, double * x, double * y, size_t length, size_t size);

int cspl_eval_periodic_min (unsigned int * min_n, double * signal, size_t size, double level);
int cspl_eval_periodic_max (unsigned int * max_n, double * signal, size_t size, double level);

/**
 * Diese Methode berechnet ein Durchschnittsfunktion aus einer Liste mit Startpunkten zur Durchschnittsbildung.
 * @author Hermann Sonntag
 * @param templ Pointer zu der Stelle, an der die Durchschnittsfunktion abgespeichert wird
 * @param signal Originalfunktion, welche die Einzelfunktionen zur Mittelung enthält
 * @param n Liste der Startpunkte
 * @param size Länge der Originalfunktion
 * @param count Anzahl der Startpunkte
 * @return Länge der gemittelten Funktion
 */
int cspl_average (double * templ, double * signal, unsigned int * n, size_t size, int count);

__END_DECLS /* Hier enden die Deklarationen für cspl. */

#endif /* __CSPL_H__ */

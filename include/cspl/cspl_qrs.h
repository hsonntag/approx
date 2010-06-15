/*   cspl/cspl_qrs.h                                                       *
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
#ifndef __CSPL_QRS_H__
#define __CSPL_QRS_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit_nlin.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS /* Hier beginnen die Deklarationen  für cspl_qrs. */

/** Beschleuniger für Interpolationen */
gsl_interp_accel * acc;

/** Parameter, die für eine Splineapproximation benötigt werden */
struct cspl_qrs_data {
    size_t n; /**< @brief Länge des Splines */
    double * t; /**< @brief Zeitvektor für die Skalierung der Zeitachse */
    double * y; /**< @brief Originalfunktion */
    gsl_spline * n_qrs; /**< @brief Spline */
    double * sigma; /**< @brief Standardabweichung */
    size_t p; /**< @brief Anzahl der Parameter, die approximiert werden sollen */
    gsl_multifit_fdfsolver * s; /**< @brief Gleichungslöser */
    gsl_matrix * covar; /**< @brief Kovarianzmatrix */
    gsl_vector_view x; /**< @brief Parametervektor mit den Initialwerten*/
};

/** Initialisierung für Splineberechnung.
 * @author Hermann Sonntag
 */
int cspl_qrs_init();

/** Deinitialisierung der Splineberechnung.
 * @author Hermann Sonntag
 */
int cspl_qrs_free();

/**
 * Berechnet Funktionswerte in Abhängigkeit des Splines und der Parameter
 * @author Hermann Sonntag
 * @param x Parametervektor
 * @param data Splineobjekt
 * @param f Funkionsvektor
 * @return Status der Funktionswerteberechnung
 */
int cspl_qrs_f (const gsl_vector * x, void * data, gsl_vector * f);

/**
 * Berechnet Jacobi-Matrix einer Funkion in Abhängigkeit des Splines und der Parameter
 * @author Hermann Sonntag
 * @param x Parametervektor
 * @param data Splineobjekt
 * @param J Jocobi-Matrix
 * @return Status der Funktionswerteberechnung
 */
int cspl_qrs_df (const gsl_vector * x, void * data, gsl_matrix * J);

/**
 * Berechnet Funktionswerte und zugehörige Jacobi-Matrix in Abhängigkeit des Splines und der Parameter
 * @author Hermann Sonntag
 * @param x Parametervektor
 * @param data Splineobjekt
 * @param f Funktionsvektor
 * @param J Jocobi-Matrix
 * @return Status der Funktionswerteberechnung
 */
int cspl_qrs_fdf (const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);

/** Enthält Funktion, Jocobi-Matrix und Parameter für eine Approximation von Parametern */
struct cspl_qrs_function_fdf_struct
{
    int (* f) (const gsl_vector * x, void * params, gsl_vector * f);
    int (* df) (const gsl_vector * x, void * params, gsl_matrix * df);
    int (* fdf) (const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix *df);
    size_t n;   /**< @brief Anzahl der Funktionswerte */
    size_t p;   /**< @brief Anzahl unabhängiger Parameter für eine Approximation*/
    void * params; /**< @brief weitere Parameter */
};

/**
 * @typedef cspl_qrs_function_fdf
 * Definition für cspl_qrs_function_fdf_struct
 */
typedef struct cspl_qrs_function_fdf_struct cspl_qrs_function_fdf ;

__END_DECLS /* Hier enden Die Deklarationen für cspl_qrs. */

#endif /* __CSPL_QRS_H__ */

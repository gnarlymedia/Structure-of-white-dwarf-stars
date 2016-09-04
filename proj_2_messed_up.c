#include <stdio.h>
#include <math.h>
#include <cpgplot.h>
#include <stdbool.h>

// global vars
#define MAX_SIZE 10000000

double Sirius_B_m_centre = 1.053;
double Sirius_B_m_error = 0.028;
double Sirius_B_r_centre = 0.0074;
double Sirius_B_r_error = 0.0006;
double Eri_B_m_centre = 0.48;
double Eri_B_m_error = 0.02;
double Eri_B_r_centre = 0.0124;
double Eri_B_r_error = 0.0005;

static double r_bar_array[MAX_SIZE];
static double rho_bar_array[MAX_SIZE];
static double m_bar_array[MAX_SIZE];

int counter;

double gamma_function(double x) {
    return x * x / (3.0 * sqrt(1.0 + pow(x, 2.0)));
}

double d_m_bar_d_r_bar(double r_bar, double rho_bar) {
    return pow(r_bar, 2.0) * rho_bar;
}

double d_rho_bar_d_r_bar(double r_bar, double rho_bar, double m_bar) {
    if (rho_bar < 0.0) {
        rho_bar = fabs(rho_bar);
    }
    return (-1.0 * m_bar * rho_bar) / (gamma_function(pow(rho_bar, 1.0 / 3.0)) * pow(r_bar, 2.0));
}

double m_bar(double h, double rho_sub_c) {
    return pow(h, 3.0) * rho_sub_c / 3.0;
}

double rho_bar(double h, double rho_sub_c) {
    return rho_sub_c * (1.0 - (pow(h, 2.0) * rho_sub_c) / (6.0 * gamma_function(pow(rho_sub_c, 1.0 / 3.0))));
}

void euler(double h, double rho_bar_sub_c, float r_bar_plot[], float rho_bar_plot[], float m_bar_plot[], double *m_bar_final,
      double *r_bar_final) {
    r_bar_array[0] = h;
    rho_bar_array[0] = rho_bar(h, rho_bar_sub_c);
    m_bar_array[0] = m_bar(h, rho_bar_sub_c);

    r_bar_plot[0] = (float) r_bar_array[0];
    rho_bar_plot[0] = (float) rho_bar_array[0];
    m_bar_plot[0] = (float) m_bar_array[0];

//    int counter = 0;
    double old_mass, mass, radius, density;

    //printf("Euler\n");
    //printf("Final values\n");
    //printf("\tRadius\t\tDensity\t\tMass\n");

    for (counter = 1; rho_bar_array[counter - 1] > 0 && counter < MAX_SIZE; counter++) {
        mass = m_bar_array[counter - 1];
        radius = r_bar_array[counter - 1];
        density = rho_bar_array[counter - 1];

//        if (counter % 10 == 0) {
//            printf("\t%lf\t%lf\t%lf\n", radius, density, mass);
//        }

        old_mass = mass;
        mass = mass + h * d_m_bar_d_r_bar(radius, density);
        density = density + h * d_rho_bar_d_r_bar(radius, density, old_mass);

        r_bar_array[counter] = radius + h;
        rho_bar_array[counter] = density;
        m_bar_array[counter] = mass;

        r_bar_plot[counter] = (float) (radius) + (float) (h);
        rho_bar_plot[counter] = (float) density;
        m_bar_plot[counter] = (float) mass;
    }

//    printf("\t%lf\t%lf\t%lf\n\n", r_bar_array[counter - 2], rho_bar_array[counter - 2], m_bar_array[counter - 2]);

    *m_bar_final = m_bar_array[counter - 1];
    *r_bar_final = r_bar_array[counter - 1];

    return;
}

double f(double x, double y, double z) {
    return d_rho_bar_d_r_bar(x, y, z);
}

double g(double x, double y, double z) {
    return d_m_bar_d_r_bar(x, y);
}

void runge_kutta(double h, double this_rho_bar_sub_c, float r_bar_plot[], float rho_bar_plot[], float m_bar_plot[],
                 double *m_bar_final, double *r_bar_final) {
    r_bar_array[0] = h;
    rho_bar_array[0] = rho_bar(h, this_rho_bar_sub_c);
    m_bar_array[0] = m_bar(h, this_rho_bar_sub_c);

    r_bar_plot[0] = (float) r_bar_array[0];
    rho_bar_plot[0] = (float) rho_bar_array[0];
    m_bar_plot[0] = (float) m_bar_array[0];

    int counter;
    double x, y, z, x_next, y_next, z_next, f_1, g_1, f_2, g_2, f_3, g_3, f_4, g_4;

    //printf("Runge-Kutta\n");
    //printf("Final values\n");
    //printf("\tRadius\t\tDensity\t\tMass\n");

    for (counter = 1; rho_bar_array[counter - 1] > 0 && counter < MAX_SIZE; counter++) {
        x = r_bar_array[counter - 1];
        y = rho_bar_array[counter - 1];
        z = m_bar_array[counter - 1];

        f_1 = f(x, y, z);
        g_1 = g(x, y, z);

        f_2 = f(
                x + h / 2.0,
                y + (h / 2.0) * f_1,
                z + (h / 2.0) * g_1
        );

        g_2 = g(
                x + h / 2.0,
                y + (h / 2.0) * f_1,
                z + (h / 2.0) * g_1
        );

        f_3 = f(
                x + h / 2.0,
                y + (h / 2.0) * f_2,
                z + (h / 2.0) * g_2
        );

        g_3 = g(
                x + h / 2.0,
                y + (h / 2.0) * f_2,
                z + (h / 2.0) * g_2
        );

        f_4 = f(
                x + h,
                y + h * f_3,
                z + h * g_3
        );

        g_4 = g(
                x + h,
                y + h * f_3,
                z + h * g_3
        );

        x_next = x + h;
        y_next = y + (h / 6.0) * (f_1 + 2.0 * f_2 + 2.0 * f_3 + f_4);
        z_next = z + (h / 6.0) * (g_1 + 2.0 * g_2 + 2.0 * g_3 + g_4);

        r_bar_array[counter] = x_next;
        rho_bar_array[counter] = y_next;
        m_bar_array[counter] = z_next;

        r_bar_plot[counter] = (float) r_bar_array[counter];
        rho_bar_plot[counter] = (float) rho_bar_array[counter];
        m_bar_plot[counter] = (float) m_bar_array[counter];
    }

    *m_bar_final = m_bar_array[counter - 1];
    *r_bar_final = r_bar_array[counter - 1];

//    printf("\t%lf\t%lf\t%lf\n", x, y, z);

    return;
}

int
plot(char *x_axis_label, char *y_axis_label, char *plot_label, double x_min, double x_max, double y_min, double y_max,
     int number_points, float x_array[], float y_array[]) {
    // printf("y_array[number_points - 1] %f", y_array[number_points - 1]);

    // set up environment
    // cpgenv(x_min, x_array[number_points - 1], y_min, y_array[number_points - 2], 0, 2);
    cpgenv(x_min, x_max, y_min, y_max, 0, 2);

    cpgbbuf();

    // set colour to black
    cpgsci(1);

    // labels
    cpglab(x_axis_label, y_axis_label, plot_label);

    cpgsci(3);

    cpgline((float) number_points, x_array, y_array);
    cpgsci(1);
    cpgebuf();
    cpgsave();

    return 0;
}

double val_from_scaled_val(double var_naught, double var_bar) {
    return var_naught * var_bar;
}

double rho_naught_func(double Y_e) {
    return 9.79E5 * pow(Y_e, -1.0);
}

double R_naught_func(double Y_e) {
    return 7.72E8 * Y_e;
}

double M_naught_func(double Y_e) {
    return 5.67E33 * pow(Y_e, 2.0);
}

double rho_solar_func(double rho) {
    return rho / 150;
}

double r_solar_func(double r) {
    return r / 6.95E10;
}

double m_solar_func(double m) {
    return m / 1.98E33;
}

int check_within_range(double lower, double upper, double val, int * results_success, char * results_success_str, size_t buff_len, char * dwarf_name) {
    * results_success = (lower < val) && (val < upper);
    if (results_success)
    {
        strncpy(results_success_str, dwarf_name, buff_len);
        results_success_str[buff_len - 1] = '\0';
    } else {}
}

int main(void) {
    printf("Start\n\n");

    double h;
    double final_m_bar_euler;
    double final_r_bar_euler;
    double final_m_bar_rk;
    double final_r_bar_rk;
    double converge_to_r_bar_rk = 1.591630;
    double converge_to_m_bar_rk = 1.298013;

    static float r_bar_array_euler_plot[MAX_SIZE];
    static float rho_bar_array_euler_plot[MAX_SIZE];
    static float m_bar_array_euler_plot[MAX_SIZE];

    static float r_bar_array_rk_plot[MAX_SIZE];
    static float rho_bar_array_rk_plot[MAX_SIZE];
    static float m_bar_array_rk_plot[MAX_SIZE];

    double m_diff_euler_pc, r_diff_euler_pc;
    double m_diff_rk_pc, r_diff_rk_pc;

    printf("Stability of solutions against step size, h for scaled central density rho_bar_sub_c = 10.0\n\n");
    printf("\t\tEuler\t\t\t\t\t\tRunge-Kutta\n");
    printf("Step\t\tM_bar\t\tR_bar\t\tM_Diff\t\tR_Diff\t\tM_bar\t\tR_bar\t\tM_Diff\t\tR_Diff\n");
    double rho_bar_sub_c = 10.0;
    for (h = 1.0; h > 1.0E-6; h = h / 10.0) {
        euler(h, rho_bar_sub_c, r_bar_array_euler_plot, rho_bar_array_euler_plot, m_bar_array_euler_plot,
              &final_m_bar_euler, &final_r_bar_euler);

        runge_kutta(h, rho_bar_sub_c, r_bar_array_rk_plot, rho_bar_array_rk_plot, m_bar_array_rk_plot, &final_m_bar_rk,
                    &final_r_bar_rk);

        r_diff_euler_pc = (final_r_bar_euler - converge_to_r_bar_rk) * 100 / converge_to_r_bar_rk;
        m_diff_euler_pc = (final_m_bar_euler - converge_to_m_bar_rk) * 100 / converge_to_m_bar_rk;

        r_diff_rk_pc = (final_r_bar_rk - converge_to_r_bar_rk) * 100 / converge_to_r_bar_rk;
        m_diff_rk_pc = (final_m_bar_rk - converge_to_m_bar_rk) * 100 / converge_to_m_bar_rk;

        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", h, final_m_bar_euler, final_r_bar_euler,
               m_diff_euler_pc, r_diff_euler_pc, final_m_bar_rk, final_r_bar_rk, m_diff_rk_pc, r_diff_rk_pc);
    }

    double rho_sub_c;

    printf("\n\nSolutions for different scaled central densities, rho_c_bar, for step size h = 0.00001 using Runge_Kutta\n\n");
    printf("rho_c_bar\t\trho_c\t\t\tM_bar\t\tR_bar\n");
    h = 0.00001;
    double Y_e = 1;
    for (rho_bar_sub_c = 0.1; rho_bar_sub_c < 1.0E9; rho_bar_sub_c = rho_bar_sub_c * 10.0) {
        runge_kutta(h, rho_bar_sub_c, r_bar_array_rk_plot, rho_bar_array_rk_plot, m_bar_array_rk_plot, &final_m_bar_rk,
                    &final_r_bar_rk);

        rho_sub_c = val_from_scaled_val(rho_naught_func(Y_e), rho_bar_sub_c);

        printf("%lf", rho_bar_sub_c);

        if (rho_bar_sub_c < 1.0E8) {
            printf("\t");
        }

        printf("\t%lf", rho_sub_c);

        if (rho_bar_sub_c < 1.0E3) {
            printf("\t");
        }

        printf("\t%lf\t%lf\n", final_m_bar_rk, final_r_bar_rk);
    }

    double r, m, rho_sub_c_solar, r_solar, m_solar;
    double Sirius_B_m_lower = Sirius_B_m_centre - Sirius_B_m_error;
    double Sirius_B_m_upper = Sirius_B_m_centre + Sirius_B_m_error;
    double Sirius_B_r_lower = Sirius_B_r_centre - Sirius_B_r_error;
    double Sirius_B_r_upper = Sirius_B_r_centre + Sirius_B_r_error;
    double Eri_B_m_lower = Eri_B_m_centre - Eri_B_m_error;
    double Eri_B_m_upper = Eri_B_m_centre + Eri_B_m_error;
    double Eri_B_r_lower = Eri_B_r_centre - Eri_B_r_error;
    double Eri_B_r_upper = Eri_B_r_centre + Eri_B_r_error;

    int results_success_Sirius_1;
    int results_success_Sirius_2;
    int results_success_Eri_1;
    int results_success_Eri_2;
    char * results_str;
    char * Sirius_m_success_str;
    char * Sirius_r_success_str;
    char * Eri_m_success_str;
    char * Eri_r_success_str;
    char * success_str;

    printf("\n\nSolutions for different scaled central densities, rho_c_bar, and Y_e, for step size h = 0.00001 using Runge_Kutta\n\n");
    printf("rho_bar_sub_c\tY_e\t\tRho_sub_c_solar\t\tM_solar\t\tR_solar\tSirius_B_Mass?\tSirius_B_Radius?\tEri_B_Mass?\tEri_B_Radius?\n");
    h = 0.0001;
    for (rho_bar_sub_c = 0.01; rho_bar_sub_c < 1.5; rho_bar_sub_c = rho_bar_sub_c + 0.01) {
        for (Y_e = 0.464; Y_e < 0.48; Y_e = Y_e + 0.001) {
            runge_kutta(h, rho_bar_sub_c, r_bar_array_rk_plot, rho_bar_array_rk_plot, m_bar_array_rk_plot,
                        &final_m_bar_rk, &final_r_bar_rk);

            rho_sub_c = val_from_scaled_val(rho_naught_func(Y_e), rho_bar_sub_c);
            r = val_from_scaled_val(R_naught_func(Y_e), final_r_bar_rk);
            m = val_from_scaled_val(M_naught_func(Y_e), final_m_bar_rk);

            rho_sub_c_solar = rho_solar_func(rho_sub_c);
            r_solar = r_solar_func(r);
            m_solar = m_solar_func(m);


            results_success_Sirius_1 = 0;
            results_success_Sirius_2 = 0;
            results_success_Eri_1 = 0;
            results_success_Eri_2 = 0;
            results_str = "";
            Sirius_m_success_str = "";
            Sirius_r_success_str = "";
            Eri_m_success_str = "";
            Eri_r_success_str = "";
            success_str = "";

            check_within_range(Sirius_B_m_lower, Sirius_B_m_upper, m_solar, &results_success_Sirius_1, success_str, sizeof(comp_name), "\t*SIRIUS B MASS");
            check_within_range(Sirius_B_r_lower, Sirius_B_r_upper, r_solar, &results_success_Sirius_2, success_str, sizeof(comp_name), "\t*SIRIUS B RADIUS");
            check_within_range(Eri_B_m_lower, Eri_B_m_upper, m_solar, &results_success_Eri_1, success_str, sizeof(comp_name), "\t*ERI B MASS");
            check_within_range(Eri_B_r_lower, Eri_B_r_upper, r_solar, &results_success_Eri_2, success_str, sizeof(comp_name), "\t*ERI B RADIUS");

            if (results_success_Sirius_1 ||
                results_success_Sirius_2 ||
                results_success_Eri_1 ||
                results_success_Eri_2
            ) {
                printf("%lf", rho_bar_sub_c);

                if (rho_bar_sub_c < 1.0E8) {
                    printf("\t");
                }

                printf("%lf", Y_e);

                printf("\t%lf", rho_sub_c_solar);

                if (rho_bar_sub_c < 1.0E5) {
                    printf("\t");
                }

                printf("\t%lf\t%lf", m_solar, r_solar);

                if ((Sirius_m_success_str != NULL) && (Sirius_m_success_str[0] != '\0')) {
                    printf("%s", Sirius_m_success_str);
                } else {
                    printf("\t\t");
                }

                if ((Sirius_r_success_str != NULL) && (Sirius_r_success_str[0] != '\0')) {
                    printf("%s", Sirius_r_success_str);
                } else {
                    printf("\t\t");
                }

                if ((Eri_m_success_str != NULL) && (Eri_m_success_str[0] != '\0')) {
                    printf("%s", Eri_m_success_str);
                } else {
                    printf("\t\t");
                }

                if ((Eri_r_success_str != NULL) && (Eri_r_success_str[0] != '\0')) {
                    printf("%s", Eri_r_success_str);
                } else {
                    printf("\t\t");
                }

                printf("\n");
            }

        }
    }

    // Plotting
//    if (1 != cpgbeg(0, "?", 1, 1))
//    if (1 != cpgbeg(0, "proj_1_plot.ps/VCPS", 1, 1))
/*    if (1 != cpgbeg(0, "/XWINDOW", 1, 1))
    {
        return 1;
    }
    // euler
    plot("Radius", "Mass", "Euler", 0.0, 1.55, 0.0, 1.35, counter, r_bar_array_euler_plot, m_bar_array_euler_plot);
    plot("Radius", "Density", "Euler", 0.0, 1.6, 0.0, 10.0, counter, r_bar_array_euler_plot, rho_bar_array_euler_plot);

    // Runge-Kutta
    plot("Radius", "Mass", "Range-Kutta", 0.0, 1.6, 0.0, 1.3, counter, r_bar_array_rk_plot, m_bar_array_rk_plot);
    plot("Radius", "Density", "Range-Kutta", 0.0, 1.6, 0.0, 10.0, counter, r_bar_array_rk_plot, rho_bar_array_rk_plot);
    cpgend();
*/

    printf("Finish\n\n");
}

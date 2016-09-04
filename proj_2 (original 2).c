#include <stdio.h>
#include <math.h>
#include <cpgplot.h>
#include <stdbool.h>

// global vars
#define MAX_SIZE 10000

double r_bar_array[MAX_SIZE];
double rho_bar_array[MAX_SIZE];
double m_bar_array[MAX_SIZE];

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

void euler(double h, double rho_bar_sub_c, float r_bar_plot[], float rho_bar_plot[], float m_bar_plot[]) {
    r_bar_array[0] = h;
    rho_bar_array[0] = rho_bar(h, rho_bar_sub_c);
    m_bar_array[0] = m_bar(h, rho_bar_sub_c);

    r_bar_plot[0] = (float)r_bar_array[0];
    rho_bar_plot[0] = (float)rho_bar_array[0];
    m_bar_plot[0] = (float)m_bar_array[0];

//    int counter = 0;
    double old_mass, mass, radius, density;

    printf("Euler\n");
    printf("Final values\n");
    printf("\tRadius\t\tDensity\t\tMass\n");

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

        r_bar_plot[counter] = (float)(radius) + (float)(h);
        rho_bar_plot[counter] = (float)density;
        m_bar_plot[counter] = (float)mass;
    }

    printf("\t%lf\t%lf\t%lf\n\n", r_bar_array[counter - 2], rho_bar_array[counter - 2], m_bar_array[counter - 2]);

    return;
}

double f(double x, double y, double z) {
    return d_rho_bar_d_r_bar(x, y, z);
}

double g(double x, double y, double z) {
    return d_m_bar_d_r_bar(x, y);
}

void runge_kutta(double h, double this_rho_bar_sub_c, float r_bar_plot[], float rho_bar_plot[], float m_bar_plot[])
{
    r_bar_array[0] = h;
    rho_bar_array[0] = rho_bar(h, this_rho_bar_sub_c);
    m_bar_array[0] = m_bar(h, this_rho_bar_sub_c);

    r_bar_plot[0] = (float)r_bar_array[0];
    rho_bar_plot[0] = (float)rho_bar_array[0];
    m_bar_plot[0] = (float)m_bar_array[0];

    int counter;
    double x, y, z, x_next, y_next, z_next, f_1, g_1, f_2, g_2, f_3, g_3, f_4, g_4;

    printf("Runge-Kutta\n");
    printf("Final values\n");
    printf("\tRadius\t\tDensity\t\tMass\n");

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

        r_bar_plot[counter] = (float)r_bar_array[counter];
        rho_bar_plot[counter] = (float)rho_bar_array[counter];
        m_bar_plot[counter] = (float)m_bar_array[counter];
    }

    printf("\t%lf\t%lf\t%lf\n", x, y, z);

    return;
}

int plot(char * x_axis_label, char * y_axis_label, char * plot_label, double x_min, double x_max, double y_min, double y_max, int number_points, float x_array[], float y_array[])
{
    // set up environment
    cpgenv(x_min, x_array[number_points - 1], y_min, y_array[number_points - 1], 0, 2);

    cpgbbuf();

    // set colour to black
    cpgsci(1);

    // labels
    cpglab(x_axis_label, y_axis_label, plot_label);

    cpgsci(3);

    cpgline((float)number_points, x_array, y_array);
    cpgsci(1);
    cpgebuf();
    cpgsave();

    return 0;
}

int main(void) {
    printf("Start\n\n");

    double h = 0.01;
    double rho_bar_sub_c = 10.0;

    float r_bar_array_euler_plot[MAX_SIZE];
    float rho_bar_array_euler_plot[MAX_SIZE];
    float m_bar_array_euler_plot[MAX_SIZE];
	euler(h, rho_bar_sub_c, r_bar_array_euler_plot, rho_bar_array_euler_plot, m_bar_array_euler_plot);

    float r_bar_array_rk_plot[MAX_SIZE];
    float rho_bar_array_rk_plot[MAX_SIZE];
    float m_bar_array_rk_plot[MAX_SIZE];
    runge_kutta(h, rho_bar_sub_c, r_bar_array_rk_plot, rho_bar_array_rk_plot, m_bar_array_rk_plot);

    // Plotting
//    if (1 != cpgbeg(0, "?", 1, 1))
//    if (1 != cpgbeg(0, "proj_1_plot.ps/VCPS", 1, 1))
    if (1 != cpgbeg(0, "/XWINDOW", 1, 1))
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

    printf("Finish\n\n");
}

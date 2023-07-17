// reference: StackOverflow question - qhull Library - C++ Interface (https://stackoverflow.com/a/72187622/6609908)

#include <iostream>
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullPoints.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <armadillo>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include "gnuplot-iostream.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>

#include <CGAL/minkowski_sum_2.h>

#include <list>
#include <cassert>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel KernelExactExact;
typedef KernelExactExact::Point_2                         Point_2;
typedef CGAL::Polygon_2<KernelExactExact>                 Polygon_2;

constexpr int dim_size = 2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel KernelExactInexact;
typedef CGAL::Convex_hull_traits_adapter_2<KernelExactInexact,
        CGAL::Pointer_property_map<Point_2>::type > Convex_hull_traits_2;

void plot_contour(arma::mat points) {
    static Gnuplot gp;

    gp << "plot '-' with linespoints pointtype 5 pointsize 3, "
          "'-' with linespoints pointtype 5 pointsize 3\n";
    gp.send1d(points);
    gp.send1d(arma::mat(points.rows(arma::uvec{{0,points.n_rows-1}})));
}

void plot_contour(const std::vector<Point_2>& points) {
    static Gnuplot gp(stdout);

    gp << "plot '-' with linespoints pointtype 5 pointsize 3" << std::endl;
    std::for_each(points.begin(), points.end(), [&](const auto &item) {
        gp << item;
    });
}

template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, std::array<T, N> arr) {

    os << "[" << arr[0];

    std::for_each(std::cbegin(arr) + 1, std::cend(arr), [&](auto e){os << ", " << e;});

    os << "]";

    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec) {

    os << "{" << vec[0] << "}";

    std::for_each(std::cbegin(vec)+1, std::cend(vec), [&](auto e){os << ", {" << e << "}";});

    return os;
}

class MyQhull{
public:
    typedef std::vector<Point_2> vecpoints ;
    typedef std::vector<std::size_t> point_indexes ;

    MyQhull(vecpoints&& vec):points(std::move(vec)), vertex_indices(), vertex_invalidated(true){}
    MyQhull(const vecpoints& vec):points(vec), vertex_indices(), vertex_invalidated(true){}

    const point_indexes& get_vertex_indices(){
        if(vertex_indices.empty() || vertex_invalidated) {
            point_indexes indices(points.size());
            std::iota(std::begin(indices), std::end(indices), 0);
            CGAL::convex_hull_2(std::begin(indices), std::end(indices),
                                std::back_inserter(vertex_indices),
                                Convex_hull_traits_2(CGAL::make_property_map(points)));
        }
        return vertex_indices;
    }

    void push_back(const Point_2& p){
        points.push_back(p);
        vertex_invalidated = true;
    }

private:
    vecpoints points;
    point_indexes vertex_indices;
    bool vertex_invalidated;
};

arma::mat armamat_from_dimension_limits(const std::vector<Point_2>& limits) {
    return arma::mat({{limits[0][0], limits[1][0]},
             {limits[0][0], limits[1][1]},
             {limits[0][1], limits[1][0]},
             {limits[0][1], limits[1][1]}});
}

double xablau(const std::vector<Point_2>& limits) {
    auto a = limits[0].x()*2.0;
    double x;
    return a.get_relative_precision_of_to_double();
}

std::vector<Point_2> vecpoints_from_dimension_limits(const std::vector<Point_2>& limits) {
    return {{{limits[0][0], limits[1][0]},
             {limits[0][0], limits[1][1]},
             {limits[0][1], limits[1][0]},
             {limits[0][1], limits[1][1]}}};
}

arma::mat vecpoint2_to_armamat(const std::vector<Point_2>& vec) {
    arma::mat m(vec.size(), 2);

    for(int i = 0; i < vec.size(); i++){
        m(i, 0) = vec[i][0];
        m(i, 1) = vec[i][1];
    }

    return m;
}

auto mRPI(const arma::mat& Acl, const std::vector<Point_2>& w, const std::vector<Point_2>& I, int max_iteration = 10) {

    std::vector<Point_2> EW = vecpoints_from_dimension_limits(w);
    Polygon_2 EWpolygon(EW.begin(), EW.end());

    std::vector<Point_2> PHIvec = vecpoints_from_dimension_limits(I);
    Polygon_2 PHIpolygon;

    std::vector<std::size_t> idx(PHIvec.size()), qhull_idx;
    std::iota(idx.begin(), idx.end(), 0);

    auto end = CGAL::convex_hull_2(idx.begin(), idx.end(), std::back_inserter(qhull_idx),
                        Convex_hull_traits_2(CGAL::make_property_map(PHIvec)));

    arma::mat PHImat(2, std::distance(qhull_idx.begin(), qhull_idx.end()));

    for (auto i = 0; i < qhull_idx.size(); i++) {
        PHImat(0, i) = PHIvec[qhull_idx[i]][0];
        PHImat(1, i) = PHIvec[qhull_idx[i]][1];
    }

    std::cout << "PHIvec:\n";
    std::cout << PHIvec << "\n";

    PHImat.print("PHImat:\n");

    plot_contour(PHImat.t());

    for(int i=0; i < max_iteration; i++) {
        PHImat = Acl * PHImat;
        PHIpolygon.
        PHImat.print("PHImat:\n");
        plot_contour(PHImat.t());
    }

    return 2;
}

int main() {

    std::cout << "mRPI program" << std::endl << std::endl << std::endl;

    // Dynamic controlled matrix (eig must lie inside the unit circle)
    arma::mat Acl{{0.879, 0.393}, {-.374, .209}};

    // max an min disturbance for each dimension: {{x1_max, x1_min}, {x2_max, x2_min}}
    std::vector<Point_2> w{{.211, -.211}, {.65, -.65}};

    // First RPI guess set vertices
    std::vector<Point_2> I = {{{5, -5}, {5, -5}}};

    // Just to see if the selected point lies inside the calculated mRPI
    Point_2 test_point{0.5,1.075};

    auto phi = mRPI(Acl, w, I, 50);

    std::cout << Acl.size() << std::endl;

    return 0;
}

#undef Qhull
#undef runQhull

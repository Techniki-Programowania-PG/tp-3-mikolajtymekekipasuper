#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "biblioteka.hpp"
#include "point.hpp"

namespace py = pybind11;

PYBIND11_MODULE(biblioteka, m)
{
      m.doc() = "Przykładowy moduł C++ z pybind11";

      py::class_<point>(m, "point")
          .def(py::init<>())
          .def(py::init<double, double>())
          .def_readwrite("x", &point::x)
          .def_readwrite("y", &point::y)
          .def("__repr__", [](const point &p)
               { return "<Point x=" + std::to_string(p.x) + " y=" + std::to_string(p.y) + ">"; });

      m.def("sinus", &sinus, "Wyznacza wartosci funkcji sinus",
            py::arg("Liczba próbek"), py::arg("Czas początkowy"), py::arg("Czas końcowy"), py::arg("Częstotliwość"));

      m.def("cosinus", &cosinus, "Wyznacza wartosci funkcji cosinus",
            py::arg("Liczba próbek"), py::arg("Czas początkowy"), py::arg("Czas końcowy"), py::arg("Częstotliwość"));

      m.def("rectangle", &rectangle, "Wyznacza wartosci sygnału prostokątnego",
            py::arg("Liczba próbek"), py::arg("Czas początkowy"), py::arg("Czas końcowy"), py::arg("Częstotliwość"));

      m.def("saw", &saw, "Wyznacza wartosci sygnału piłokształtnego",
            py::arg("Liczba próbek"), py::arg("Czas początkowy"), py::arg("Czas końcowy"), py::arg("Częstotliwość"));

      m.def("draw", &draw, "Rysuje wykres funkcji",
            py::arg("Dane"), py::arg("Opis wykresu"));

      m.def("DFT", &DFT, "Dyskretna Transformata Fouriera",
            py::arg("Sygnał wejściowy"));
      m.def("spectrumDFT", &spectrumDFT, "Widmo DFT",
            py::arg("DFT sygnału"));
      m.def("IDFT", &IDFT, "Odwrotna Dyskretna Transformata Fouriera",
            py::arg("DFT sygnału"), py::arg("Czas początkowy"), py::arg("Czas końcowy"));
      m.def("filter1D", &filter1D, "Filtr 1D",
            py::arg("DFT sygnału"), py::arg("Modyfikator"));
      m.def("generateMatrix", &generateMatrix, "Generator Macierzy",
            py::arg("Wiersze"), py::arg("Kolumny"));
      m.def("applyGaussianBlur2D", &applyGaussianBlur2D, "Rozmycie Gaussa",
            py::arg("Sygnał wejściowy"), py::arg("Rozmiar"), py::arg("Sigma"), py::arg("Skala"));
      m.def("plotMatrix", &plotMatrix, "Rysuj Macierz",
            py::arg("Macierz"), py::arg("Opis wykresu"));
      m.def("convert1Dto2D", &convert1Dto2D, "1D do 2D",
            py::arg("Sygnał wejściowy"), py::arg("Skala"));
}
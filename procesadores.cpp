/*
* procesadores.cpp
*
*  Created on: 21/4/2015
*      Author: Miguel Ascanio Gómez
*      Universidad Complutense de Madrid
*/

// Compilación gcc (versiones recientes)
// g++ -O3 -std=c++11 -o procesadores procesadores.cpp
// Compilación gcc (versiones anteriores)
// g++ -O3 -std=c++0x -o procesadores procesadores.cpp
// Compilación MinGW (existe un bug en math.h)
// g++ -O3 -std=c++0x -D__NO_INLINE__ -o procesadores procesadores.cpp

// Comentar para usar la implementación de la STL de heap y vector, en lugar de la de BOOST
// (lo que viene a ser el montículo de fibonacci)
// #define USE_BOOST

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>

#ifdef USE_BOOST
#ifdef _WIN32
#include <boost\heap\fibonacci_heap.hpp>
#include <boost\container\vector.hpp>
#else
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/container/vector.hpp>
#endif
#else
#include <queue>
#include <vector>
#endif

#include <algorithm>
#include <functional>

#include <ctime>
#include <stdint.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#endif

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// VARIABLES GLOBALES /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

uint32_t NumTask = 15;	// Num tareas (valor por defecto)
uint32_t NumProcs = 4;	// Num procesadores (valor por defecto)

uint32_t iters = 20;	// Numero de ejecuciones de cada test (para genera test)  (valor por defecto)

bool silent = true;		// Para quitar parte de la salida. Algunas opciones de ejecución modifican el valor.

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// CONSTANTES /////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

const uint32_t MaxNumProcs = 8;		// Máximo número procesadores

const uint32_t SEED = 1;	

const uint32_t T_TASK_MIN = 1;
const uint32_t T_TASK_MAX = 10;
const uint32_t T_TASK_STEP = 2;

#define USHRT_MAX 0xffff

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// DEFINICIONES DE TIPOS //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

typedef uint8_t procId_t;	// Tipo para el identificador del procesador
typedef uint8_t taskId_t;	// Tipo para el identificador de tarea
typedef uint16_t taskt_t;	// Tipo para el tiempo de las tareas (tanto acumulado como de cada tarea individualmente)

const taskt_t INF = 0xffff; // USHORT MAX, debería cambiarse para estar acorde al tipo de taskt_t

#ifdef USE_BOOST
typedef boost::container::vector <taskId_t> vtaskId_t;
typedef boost::container::vector <taskt_t> vtaskt_t;
#else
typedef std::vector<taskId_t> vtaskId_t;
typedef std::vector<taskt_t> vtaskt_t;
#endif
// Estructura que representa a un procesador, esto es, la carga que tiene de tareas y su id
struct processor_t {
	taskt_t tiempoAcuMax;
	procId_t id;

	processor_t() {}

	processor_t (procId_t a) {
		id = a;
		tiempoAcuMax = 0;
	}
};

bool operator< (const processor_t& lhs, const processor_t& rhs) {
	if ( lhs.tiempoAcuMax != rhs.tiempoAcuMax)
		return lhs.tiempoAcuMax < rhs.tiempoAcuMax;
	else
		return lhs.id < rhs.id;
}

bool operator> (const processor_t& lhs, const processor_t& rhs) {
	if ( lhs.tiempoAcuMax != rhs.tiempoAcuMax)
		return lhs.tiempoAcuMax > rhs.tiempoAcuMax;
	else
		return lhs.id > rhs.id;
}

// Clase que representa al conjunto de los procesadores, esto es, un array de NumProcs 
// de procesadores, además de algún dato auxiliar
class processors_t {
private:
	processor_t _v[MaxNumProcs];
	procId_t posOld;
	taskt_t tOld;

public:
	processors_t() {
		for (procId_t i = 0; i < NumProcs; i++)
			_v[i] = processor_t(i);
	}

	// Copy assignment
	processors_t(const processors_t& ps) {
		std::copy(ps.begin(), ps.end(), this->_v);
		this->posOld = ps.posOld;
		this->tOld = ps.tOld;
	}

	typedef const processor_t* iterator;

	// Devuelve el procesador con id i
	const processor_t& operator[] (procId_t x) const { 
		for(int i = 0; i < NumProcs; i++) {
			if (_v[i].id == x)
				return _v[i];
		}
	}

	const processor_t* begin() const { return (std::begin(_v)); }
	const processor_t* end() const { return (_v + NumProcs); }

	// Esta funcion incrementa el valor del tiempo acumulado del procesador n en una cantidad de tiempo t
	// Devuelve el elemento cambiado y matiene el vector ordenado, MAX A LA DERECHA
	// Mayores a la derecha
	void incrementar(procId_t n, taskt_t t) {
		int i;
		for(i = 0; i < NumProcs; i++) {
			if (_v[i].id == n)
				break;
		}

		tOld = _v[i].tiempoAcuMax;
		_v[i].tiempoAcuMax += t;
		posOld = i; // Por si es la última
		// Llevar a la derecha este hijo
		for (procId_t j = i; j < NumProcs-1; j++)
		{
			if (_v[j] < _v[j+1])
				break;

			std::swap(_v[j], _v[j+1]);
			posOld = j+1;
		}
	}

	// Utilizado para la voraz optimista
	void incrementarIni(taskt_t t) {
		tOld = _v[0].tiempoAcuMax;
		_v[0].tiempoAcuMax += t;

		// Llevar a la derecha este hijo
		for (procId_t j = 0; j < NumProcs-1; j++)
		{
			if (_v[j] < _v[j+1])
				break;

			std::swap(_v[j], _v[j+1]);
			posOld = j+1;
		}
	}

	// Devuelve el último procesador modificado a su estado anterior
	// Usado tanto en la voraz optimista como al acabar de procesar un nodo
	// para generar los hijos correctamente
	void decrementar() {
		_v[posOld].tiempoAcuMax = tOld;
		// Llevar a la izquierda
		for (procId_t i = posOld; i > 0; i--) {
			if (_v[i-1] < _v[i])
				break;
			std::swap(_v[i], _v[i-1]);
		}
	}
};

// Estructura que representa a un nodo. Se guardan la solución hasta el momento, 
// los tiempos de cada procesador, la profuncidad, el coste actual y la estimación
struct node_t {
	vtaskId_t sol;
	processors_t tiemposProcs;

	taskId_t k;				// Profundidad, en la práctica (número de tareas - 1) que hemos metido
	taskt_t tiempoAcuMax;	// Coste, la carga del procesador más cargado.
	taskt_t estimacion;		// Coste optimista (prioridad)

	node_t() {
		tiemposProcs = processors_t();
		sol = vtaskId_t(NumTask, -1);
		k = -1;
		tiempoAcuMax = 0;
		estimacion = 0;
	}

	// Copy assignment
	node_t& operator= (const node_t& rhs) {
		// Copiar los vectores (el contenido)
		this->sol = rhs.sol;
		this->tiemposProcs = rhs.tiemposProcs; 

		this->k = rhs.k;
		this->tiempoAcuMax = rhs.tiempoAcuMax;
		this->estimacion = rhs.estimacion;

		return *this;
	}
};

inline bool operator> (const node_t& lhs, const node_t& rhs) {
	return (lhs.estimacion > rhs.estimacion);
}

inline bool operator< (const node_t& lhs, const node_t& rhs) {
	return (lhs.estimacion < rhs.estimacion);
}

#ifdef USE_BOOST
typedef boost::heap::fibonacci_heap<node_t,boost::heap::compare<std::greater<node_t>>> cola_t;
#else
typedef std::priority_queue<node_t, std::vector<node_t>, std::greater<node_t>> cola_t;
#endif

vtaskt_t T;		 // Variable global para los tiempos de las tareas
vtaskt_t T_REST; // Variable global para el tiempo de proceso restante (T_REST(0) = SUM(T), T_REST(end) = T[0]), size(T) + 1 = size(T_REST)

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// FUNCIONES DE ESTIMACION OPTIMISTA //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Función lambda que toma un node_t y devuelve un taskt_t (la estimación)
typedef std::function<taskt_t (const node_t&)> f_opt_t; 

// Función que no estima nada, siempre devuelve 0. Es como no usar cota optimista.
f_opt_t opt_muy_ingenua = [](const node_t& n) -> taskt_t { return 0; };
// Función muy barata que estima poco, pero bastante más que la de antes.
f_opt_t opt_ingenua = [](const node_t& n) -> taskt_t { return (*(n.tiemposProcs.end() - 1)).tiempoAcuMax; };

// Estumación similar a la voraz, iguala todos los procesadores, y reparte de manera uniforme el trabajo restante
// (la idea es que al final de la ejecución todos los procesadores tengan la misma carga, pero partiendo de la carga 
// de los procesadores del nodo n, y no desde el inicio como en la función anterior). Si no hay tareas suficientes
// para equilibrar la carga, i.e. el procesador con más carga difiere del siguiente con más carga más que 
// la suma de los tiempos de las atreas restantes, se devuelve el tiempo del procesador con más carga
f_opt_t opt_voraz = [](const node_t& n) -> taskt_t { 
	taskt_t dAcu = 0;
	taskt_t d;

	procId_t i = 1;
	for (processors_t::iterator it = n.tiemposProcs.begin(); it != (n.tiemposProcs.end()-1); it++) {
		d = (*(it+1)).tiempoAcuMax - (*it).tiempoAcuMax; // Distancia entre it y el anterior
		d = d * i; 	// Distancia entre el anterior a it y todos los siguientes
		dAcu += d;	// Distancias acumuladas		

		i++;
	}

	taskt_t tiempoTareas = T_REST[n.k+1];
	taskt_t r;
	if (tiempoTareas > dAcu)
		r = ((tiempoTareas - dAcu) / NumProcs) + (*(n.tiemposProcs.end() - 1)).tiempoAcuMax;
	else
		r = (*(n.tiemposProcs.end() - 1)).tiempoAcuMax;

	return r;
	// Necesariamente, si no me da con las tareas que me quedan a rellenar el espacio entre 
	// el procesador más cargado y los siguientes, no voy a tardar más que el procesador con máxima carga
};

taskId_t umbral_Din_Opt = 14;

f_opt_t opt_dinamica = [](const node_t& n) -> taskt_t { 
	if (n.k < umbral_Din_Opt) {
		return opt_voraz(n);
	} else {
		return opt_muy_ingenua(n);
	}
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// FUNCIONES DE ESTIMACION PESIMISTA //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Función lambda que toma un node_t y devuelve un taskt_t (la estimación)
typedef std::function<taskt_t (const node_t&, taskt_t tiempoMejor)> f_pes_t; 

// Función que no estima nada, devuelve tiempoMejor. Es como no usar cota pesimista.
f_pes_t pes_ingenua = [](const node_t& n, taskt_t tiempoMejor) -> taskt_t { return tiempoMejor; };

// Función que va colocando las tareas más grandes en los procesadores menos cargados.
// Pre: T ordenado creciente
f_pes_t pes_voraz = [](const node_t& n, taskt_t tiempoMejor) -> taskt_t {

	if ((*(n.tiemposProcs.end()-1)).tiempoAcuMax >= tiempoMejor) // Parar si acabo superando el actual
		return tiempoMejor;

	processors_t procs = processors_t(n.tiemposProcs);

	for (int i = NumTask - 1; i >= n.k; i--) {
		procs.incrementarIni(T[i]);
		if ((*(procs.end()-1)).tiempoAcuMax > tiempoMejor) // Parar si acabo superando el actual
			return tiempoMejor;
	}

	return (*(procs.end()-1)).tiempoAcuMax;
};
taskId_t umbral_Din_Pes = 2;

f_pes_t pes_dinamica = [](const node_t& n, taskt_t tiempoMejor) -> taskt_t { 
	if (n.k < umbral_Din_Pes) {
		return pes_voraz(n, tiempoMejor);
	} else {
		return pes_ingenua(n, tiempoMejor);
	}
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// FUNCIONES DE NO-REPETICION /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Función lambda que toma un nodo y la profundidad k, y decide si hay repeticiones a partir del nodo
// Sólo con el nodo y la id valdría, pero con los otros parámetros es un poco más eficiente
typedef std::function<bool (const node_t&, procId_t id, procId_t i, taskt_t cargaLastProcesIn)> f_noRep_t;

// Función que no poda nada (permite cualquier repeticion)
f_noRep_t noR_nula = [](const node_t&, procId_t, procId_t, taskt_t) -> bool { return true; };

// Función que evita las repeticiones causadas por las M primeras tareas
// Sólo permito ejecutar si la id del procesador es menor que la profundidad, 
// i.e. en profundidad 3, sólo puedo meter la tarea en los 3 primeros procesadores (o si hay menos de 3 procesadores, en todos)
f_noRep_t noR_ini  = [](const node_t& X, procId_t id, procId_t, taskt_t) -> bool { return id <= X.k; };

// Función que evita las repeticiones causadas por procesadores con la misma carga
// Sólo permito meter una tarea, en el procesador que ocupe la posición i de los procesadores de X 
// (lo que viene a ser el siguiente procesador que estoy probando, el que menos carga tiene), 
// si este procesador i tiene diferente carga que el último procesador que recibió una tarea
// (la carga es antes de recibir la tarea), esta carga es pasada como parámetro
f_noRep_t noR_mism = [](const node_t& X, procId_t id, procId_t i, taskt_t cargaLastProcesIn) -> bool { 
	return X.tiemposProcs[i].tiempoAcuMax != cargaLastProcesIn;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Tipos para los tests
#ifdef USE_BOOST
typedef boost::container::vector<f_opt_t> v_f_opt_t;
typedef boost::container::vector<f_pes_t> v_f_pes_t;
typedef boost::container::vector<f_noRep_t> v_noRep_t;
typedef boost::container::vector<uint32_t> v_seed_t;
#else
typedef std::vector<f_opt_t> v_f_opt_t;
typedef std::vector<f_pes_t> v_f_pes_t;
typedef std::vector<f_noRep_t> v_noRep_t;
typedef std::vector<uint32_t> v_seed_t;
#endif

void tareas_rp(vtaskId_t& sol_mejor, taskt_t& tiempo_mejor, uint64_t& expandedNodes, uint64_t& maxExpandedNodes, 
			   f_opt_t estimacionOptimista, f_pes_t estimacionPesimista, f_noRep_t funcionNoRep);

void generaTareas(uint32_t seed);
void generaTest(v_f_opt_t opts, v_f_pes_t pess, v_noRep_t reps, int iters, std::ostream* out, bool draw);

uint64_t test(f_opt_t estimacionOptimista, f_pes_t estimacionPesimista, f_noRep_t funcionNoRep);
bool testVal(f_opt_t estimacionOptimista, f_pes_t estimacionPesimista, f_noRep_t funcionNoRep);
void testDin();

void draw(const vtaskId_t& sol);

int ejecucionPorDefecto(v_f_opt_t opts, v_f_pes_t pess, v_noRep_t reps);
int ejecucionParametrica(v_f_opt_t opts, v_f_pes_t pess, v_noRep_t reps, 
						 uint32_t StartNumTask, uint32_t EndNumTask, uint32_t StartNumProcs, uint32_t EndNumProcs, int tipoSalida);

uint64_t GetTimeMs64();

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	int r;

	uint32_t StartNumTask;	// Num tareas iniciales
	uint32_t StartNumProcs;	// Num procesadores iniciales

	uint32_t EndNumTask;	// Num tareas finales
	uint32_t EndNumProcs;	// Num procesadores finales

	// Test a realizar
	v_f_opt_t opts;			// Funciones de estiamción optimista a probar
	v_f_pes_t pess;			// Funciones de estiamción pesimista a probar
	v_noRep_t reps;			// Funciones de no repetición a probar

	// Introducimos en los vectores los test que queremos pasar
	// Los vectores de funciones deben ser del mismo tamaño, pues se realizará una llamada
	// a tareas RP con N test introducidos en cada vector, de tal forma
	// que la i-ésima llamada a tareas RP se realizará con las funciones i-ésimas de cada vector

	opts.push_back(opt_voraz);
	pess.push_back(pes_voraz);
	reps.push_back(noR_mism);
	
	opts.push_back(opt_ingenua);
	pess.push_back(pes_voraz);
	reps.push_back(noR_mism);

	opts.push_back(opt_muy_ingenua);
	pess.push_back(pes_voraz);
	reps.push_back(noR_mism);
	
	opts.push_back(opt_voraz);
	pess.push_back(pes_ingenua);
	reps.push_back(noR_mism);

	opts.push_back(opt_ingenua);
	pess.push_back(pes_ingenua);
	reps.push_back(noR_mism);

	opts.push_back(opt_muy_ingenua);
	pess.push_back(pes_ingenua);
	reps.push_back(noR_mism);
	
	srand(SEED);

	if (argc == 1) {
		r = ejecucionPorDefecto(opts, pess, reps);
	} else if (argc == 2) {
		if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
			std::cout << "procesadores.exe" << std::endl
				<< "Usage:\tprocesadores.exe - ejecuta las opciones por defecto, salida por pantalla" << std::endl
				<< "   procesadores.exe [-s|-ss|-a|-d] iniTask endTask iniProcs endProcs iters" << std::endl
				<< "   Con: -s para salida estándar, -ss para salida estándar sólo con salida de las métricas" << std::endl
				<< "        -a para salida a archivos, -d para estandar dibujando solucion" << std::endl
				<< "   posteriormente se ejecuatran las pruebas utilizando tareas entre [iniTask, endTask]" << std::endl
				<< "   utilizando para ello [iniProcs, endProcs] para cada prueba, ejecutando cada prueba iters veces." << std::endl;
			return EXIT_FAILURE;
		}
	} else {
		if (argc != 7) {
			std::cerr << "Error, num parámetros incorrecto" << std::endl;
			return EXIT_FAILURE;
		}

		int tipoSalida;

		if (strcmp(argv[1], "-s") == 0) {
			tipoSalida = 1;
		} else if (strcmp(argv[1], "-a") == 0) {
			tipoSalida = 2;
		} else if (strcmp(argv[1], "-d") == 0) {
			tipoSalida = 3;
		} else if (strcmp(argv[1], "-ss") == 0) {
			tipoSalida = 1;
			silent = true;
		} else {
			std::cerr << "Error, el parámetro de tipo de salida no es -s o -a" << std::endl;
			return EXIT_FAILURE;
		}

		StartNumTask = atoi(argv[2]);
		EndNumTask = atoi(argv[3]);
		StartNumProcs = atoi(argv[4]);
		EndNumProcs = atoi(argv[5]);

		if (atoi(argv[4]) > MaxNumProcs) {
			std::cerr << "Error, Se excede el máximo de procesadores" << std::endl;
			return EXIT_FAILURE;
		}

		iters = atoi(argv[6]);

		r = ejecucionParametrica(opts, pess, reps, StartNumTask, EndNumTask, StartNumProcs, EndNumProcs, tipoSalida);
	}

	return r;
}

// Función que ejecuta el algoritmo de vuelta atrás con poda, utilizando las funciones pasadas como parámetro para las podas
void tareas_rp(vtaskId_t& sol_mejor, taskt_t& tiempo_mejor, uint64_t& expandedNodes, uint64_t& maxExpandedNodes, f_opt_t estimacionOptimista, f_pes_t estimacionPesimista, f_noRep_t funcionNoRep) {

	expandedNodes = 1;
	maxExpandedNodes = 1;
	uint64_t actualNodes = 1;

	cola_t cola;
	node_t X;
	taskt_t oldAcu;				// Carga del procesador antes de meterle una tarea nueva (utilizado para reestablecerlo)
	taskt_t lastProcesIn = -1;	// Carga (previa a nueva tarea) del último procesador que recibio una tarea
	taskt_t tAux; 
	
	bool tiempoMejorEsReal = false;

	cola.push(X);
	tiempo_mejor = INF; // Como infinito...

	while (!cola.empty() && (((X = cola.top()).estimacion < tiempo_mejor) || (X.estimacion == tiempo_mejor && !tiempoMejorEsReal))) {
		cola.pop();
		actualNodes--;
		X.k = X.k + 1;

		lastProcesIn = -1;	// Resetear último en entrar

		processor_t processor;
		// Probar colocando la tarea en los NumProcs procesadores
		for (procId_t i = 0; i < NumProcs; i++) 
		{
			processor = X.tiemposProcs[i];
			if (funcionNoRep(X, processor.id, i, lastProcesIn)) {
				// Cualquier camino es factible, pues la condición de ser solución es meter todas las tareas, k = NumTask
				// Pongo la tarea que estoy probando en el procesador procId
				X.sol[X.k] = processor.id;

				// Cargo el procesador con dicha tarea, guardo la carga para evitar repeticiones (si procede)
				//X.tiemposProcs[i] += T[X.k];	
				lastProcesIn = processor.tiempoAcuMax;
				X.tiemposProcs.incrementar(i, T[X.k]);

				// Guardo la carga máxima (para reestablecerla después de encolarlo)
				oldAcu = X.tiempoAcuMax;

				// Calculo el tiempo acumulado máximo y la nueva estimación
				X.estimacion = estimacionOptimista(X);
				X.tiempoAcuMax = (*(X.tiemposProcs.end() - 1)).tiempoAcuMax;

				if (X.k == NumTask-1) { // Es solución
					if ((X.tiempoAcuMax < tiempo_mejor) || (X.tiempoAcuMax == tiempo_mejor && !tiempoMejorEsReal)) {
						// Es mejor estrictamente que la última solución, o bien es igual que el último tiempo mejor
						// y además este tiempo mejor era una estimación (evitar cambiar una solución por otra de mismo coste)
						sol_mejor = X.sol;
						tiempo_mejor = X.tiempoAcuMax;
						tiempoMejorEsReal = true;
					}
				// No es solución, encolar
				} else if ((X.estimacion < tiempo_mejor) || (X.estimacion == tiempo_mejor && !tiempoMejorEsReal)) { 
					// Es mejor estrictamente que la última solución, o bien es igual que el último tiempo mejor
					// y además este tiempo mejor era una estimación (evitar cambiar una solución por otra de mismo coste)
					cola.push(X);

					tAux = estimacionPesimista(X, tiempo_mejor);
					if (tAux <= tiempo_mejor) {
						tiempo_mejor = tAux;
						tiempoMejorEsReal = false;
					}

					expandedNodes++;
					actualNodes++;
					if (actualNodes > maxExpandedNodes)
						maxExpandedNodes = actualNodes;
				}

				// Revertir cambios tras encolarlo (para el siguiente hijo)
				// X.tiemposProcs[i] -= T[X.k];
				X.tiemposProcs.decrementar();
				X.tiempoAcuMax = oldAcu;
			}
		}
	}
}

// Función que genera las tareas para el proceso, utilizando la semilla seed
void generaTareas(uint32_t seed) {
	srand(seed);
	T.clear();
	T_REST.clear();
	// Rellenar la duración de las tareas
	for (uint32_t i = 0; i < NumTask; ++i) {
		T.push_back(((rand() % (T_TASK_MAX - T_TASK_MIN)) + T_TASK_MIN) * T_TASK_STEP);
		//T.push_back(i+1);
	}
	// Ordenar las tareas (no hay perdida de generalidad)
	std::sort(T.begin(), T.end());
	// Rellenar el vector de tiempoRestante
	time_t acu = 0;
	for (uint32_t i = 0; i < NumTask; i++) {
		acu += T[i];
	}

	T_REST.push_back(acu);
	for (vtaskt_t::iterator it = T.begin(); it != T.end(); it++) {
		acu -= (*it);
		T_REST.push_back(acu);
	}

	if (!silent) {
		std::cout << std::endl << "Tiempos de ejecución de cada tarea: " << std::endl; 
		for (size_t i = 0; i < T.size(); i++) {
			std::cout << i << ": " << (int(T[i])) << std::endl ;
		}
		std::cout << std::endl << "Tiempos restantes: " << std::endl; 
		for (size_t i = 0; i < T_REST.size(); i++) {
			std::cout << i << ": " << (int(T_REST[i])) << std::endl ;
		}
	}
}

// Función que llama al algoritmo de vuelta atrás con poda, utilizando las funciones pasadas como parámetro para las podas
uint64_t test(f_opt_t estimacionOptimista, f_pes_t estimacionPesimista, f_noRep_t funcionNoRep) {
	vtaskId_t solMejor;
	taskt_t tiempoTareas = INF;
	uint64_t maxExpandedNodes, expandedNodes;

	uint64_t ini = GetTimeMs64();
	tareas_rp(solMejor, tiempoTareas, expandedNodes, maxExpandedNodes, estimacionOptimista, estimacionPesimista, funcionNoRep);
	uint64_t end = GetTimeMs64();
	if(!silent) {
		std::cout << std::endl << "Nodos expandidos: " << expandedNodes << std::endl;
		std::cout << "Máximo numero de nodos simultáneos: " << maxExpandedNodes << std::endl;
		std::cout << "Tiempo de ejecución(ms): " << end - ini << std::endl << std::endl ;
		draw(solMejor);
		std::cout << std::endl << "Tiempo maximo (tareas): " << tiempoTareas << std::endl << std::flush;
	}

	return end - ini;
}

// Función que realiza, para los elementos de los vectores de funciones optimistas, pesimistas y de no-repeticion,
// iters iteraciones con las semillas de seeds. Nótese que se llama a tareas_rp(...opts[k], pess[k], reps[k]...) iters veces, 
// no a todas las combinaciones posibles entre opts, pess y reps.
void generaTest(v_f_opt_t opts, v_f_pes_t pess, v_noRep_t reps, int iters, std::ostream* out, bool drawSol) {
	if (opts.size() != pess.size() || opts.size() != reps.size() || opts.size() == 0) {
		std::cerr << "Error, los tamaños de los vectores de prueba no coinciden o son nulos" << std::endl;
		return;
	}
	v_seed_t seeds;
	for (int i = 0; i < iters; i++) {
		seeds.push_back(rand());
	}

	vtaskId_t solMejor;
	taskt_t tiempoTareas = INF;
	uint64_t maxExpandedNodes, expandedNodes;

	for (int i = 0; i < opts.size(); i++) {
		for (int k = 0; k < iters; k++) {
			generaTareas(seeds[k]);
			uint64_t ini = GetTimeMs64();
			tareas_rp(solMejor, tiempoTareas, expandedNodes, maxExpandedNodes, opts[i], pess[i], reps[i]);
			uint64_t end = GetTimeMs64();
			(*out) << end - ini << " " << expandedNodes << " " << maxExpandedNodes << std::endl;
			if (drawSol) {
				draw(solMejor);
			}
		}
		*out << std::endl;
	}
}

// Esta función imprime por salida estándar la representación, como diagrama de Gantt, la distribución de las tareas.
void draw(const vtaskId_t& sol) {

	std::cout << "Solución: " << std::endl ; 

	for (procId_t i = 0; i < sol.size(); i++)
		std::cout << (int)sol.at(i) << " ";

	bool parity = true;
	std::cout << std::endl << "Gantt: "<< std::endl;
	for (procId_t i = 0; i < NumProcs; i++)
	{
		parity = true;
		std::cout << (int) i << ": " ;
		for (size_t j = 0; j < sol.size(); j++)
		{
			if (sol.at(j) == i) {
				for (int k = 0; k < T[j]; k++)
					std::cout << (parity ? "+" : "-");
				parity = !parity;				
			}
		}
		std::cout << std::endl;
	}
}

// Función que ejecuta los test con las opciones por defecto (llama a generaTest)
int ejecucionPorDefecto(v_f_opt_t opts, v_f_pes_t pess, v_noRep_t reps) {

	generaTest(opts, pess, reps, iters,  &(std::cout), true);

	if(!silent) {
		std::cerr << "Presione enter para continuar..." ;
		std::cin.get();
	}

	return EXIT_SUCCESS;
}

// Función que ejecuta los test con las opciones pasadas como parámetro (llama a generaTest)
int ejecucionParametrica(v_f_opt_t opts, v_f_pes_t pess, v_noRep_t reps, uint32_t StartNumTask, uint32_t EndNumTask, uint32_t StartNumProcs, uint32_t EndNumProcs, int tipoSalida) {

	for (int i = StartNumProcs; i <= EndNumProcs; i++)
	{
		NumProcs = i;
		for (int j = StartNumTask; j <= EndNumTask; j++)
		{
			NumTask = j;
			if (tipoSalida == 1) {
				generaTest(opts, pess, reps, iters, &(std::cout), false);
			} else if (tipoSalida == 2) {
				std::ofstream of;
				std::stringstream ss;
				if (i < 10)
					ss << "pr_0" << i;
				else
					ss << "pr_" << i;
				if (j < 10) 
					ss << "_task_0" << j << ".txt";
				else
					ss << "_task_" << j << ".txt";

				of.open(ss.str());
				if (!of.is_open())
					return EXIT_FAILURE;
				std::cout << "Procesadores: " << i << " Tareas: " << j << std::endl;
				generaTest(opts, pess, reps, iters, &of, false);
				of.close();
			} else {
				generaTest(opts, pess, reps, iters, &(std::cout), true);
			}
		}
	}	

	if(!silent) {
		std::cerr << "Presione enter para continuar..." ;
		std::cin.get();
	}

	return EXIT_SUCCESS;
}

// Función que compara la salida (el tiempo máximo, no la organización de als tareas) de ejecutar tareas_rp con las funciones pasadas como parámetro,
// con ejecutar tareas_rp con las funciones nulas. Algo así como un """test""" de funcionamiento.
bool testVal(f_opt_t estimacionOptimista, f_pes_t estimacionPesimista, f_noRep_t funcionNoRep) {

	generaTareas(time(NULL));

	vtaskId_t solMejor;
	taskt_t tiempoTareas = INF;
	uint64_t maxExpandedNodes, expandedNodes;

	tareas_rp(solMejor, tiempoTareas, expandedNodes, maxExpandedNodes, estimacionOptimista, estimacionPesimista, funcionNoRep);

	vtaskId_t solMejor1;
	taskt_t tiempoTareas1 = INF;

	tareas_rp(solMejor1, tiempoTareas1, expandedNodes, maxExpandedNodes, opt_muy_ingenua, pes_ingenua, noR_nula);

	draw(solMejor);
	draw(solMejor1);

	return tiempoTareas == tiempoTareas1;
}

// Función que ejecuta las versiones dinamicas de las funciones de estimación.
void testDin() {

	v_f_opt_t opts;
	v_f_pes_t pess;
	v_noRep_t reps;
	v_seed_t seeds;

	opts.push_back(opt_dinamica);
	pess.push_back(pes_voraz);
	reps.push_back(noR_mism);

	opts.push_back(opt_voraz);
	pess.push_back(pes_dinamica);
	reps.push_back(noR_mism);

	opts.push_back(opt_dinamica);
	pess.push_back(pes_dinamica);
	reps.push_back(noR_mism);

	for (int i = 0; i <= NumTask; i++) {
		umbral_Din_Opt = i;
		umbral_Din_Pes = i;

		generaTest(opts, pess, reps, iters,  &(std::cout), false);
	}
}

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
* windows and linux. */

uint64_t GetTimeMs64()
{
#ifdef _WIN32
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	* to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	uint64_t ret = li.QuadPart;
	ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
	ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

	return ret;
#else
	/* Linux */
	struct timeval tv;

	gettimeofday(&tv, NULL);

	uint64_t ret = tv.tv_usec;
	/* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
	ret /= 1000;

	/* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
	ret += (tv.tv_sec * 1000);

	return ret;
#endif
}

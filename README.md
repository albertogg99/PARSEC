# PARSEC_TFG
Versión modernizada del *benchmark* PARSEC. Este repositorio incluye 8 de las aplicaciones del *benchmark* original, las cuales han sido modernizadas durante mi Trabajo Fin de Grado: ***Estudio de la calidad de software en benchmarks científicos***. Para descargar la versión original visita https://parsec.cs.princeton.edu/

## Aplicaciones modernizadas
- Blackscholes
- Bodytrack
- Canneal
- Fluidanimate
- Freqmine
- Netstreamcluster
- Streamcluster
- Swaptions

Todas estas aplicaciones compilan sin avisos ni errores con GCC versión 10.2 y el estándar **C++17**. 

*Nota: algunas de las aplicaciones presentes en este proyecto necesitan de una librería para su funcionamiento. Dichas librerías han sido incluidas en el repositorio para permitir la ejecución del benchmark pero no han sido modernizadas, por lo que su compilación si producirá avisos (que no errores)*.

## Instalación y ejecución

A continuación se explica brevemente cómo instalar y ejecutar las aplicaciones del *benchmark* PARSEC. Para una información más detallada visita https://parsec.cs.princeton.edu/parsec3-doc.htm

#### Cómo instalar  y ejecutar el *benchmark*
**1.**  Descarga el zip con los archivos.

**2.**  Extrae los archivos en un directorio.

**3.**  Navega al directorio donde se han extraído los archivos.

**4.**  Ejecuta el comando **source env.sh** en la terminal para activar el entorno.

**5.**  Ejecuta el comando **parsecmgmt -a build** para construir las aplicaciones del benchmark.

**6.**  Ejecuta el comando **parsecmgmt -a run -i native** para ejecutar las aplicaciones del benchmark con el input nativo.

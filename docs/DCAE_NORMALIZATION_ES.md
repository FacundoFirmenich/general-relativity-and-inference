# DCAE: normalización y umbrales

**Umbrales impuestos por la normalización.** Definiendo $R=\mathcal A(s)/\mathcal A(s_0)>0$, el logaritmo se cancela en $I=\max\{0,1-R^{-1}\}$. Para $0<c<1$, $I=c$ equivale a $R=(1-c)^{-1}$ y $\Delta=-\ln(1-c)$. En particular, $I=0.90$ corresponde exactamente a $R=10$ y $\Delta=\ln10$, mientras que $I=0.95$ corresponde a $R=20$. Estos valores son consecuencias de la definición, no evidencia independiente de una barrera informacional ni de una transición física. Deben comunicarse también $R$ y $\Delta$, con referencia, incertidumbres y covarianzas; la proximidad de $I$ a uno no mide confianza ni significación estadística.

## Normalización, cambio de base y umbrales construidos

Para $R=\mathcal A(s)/\mathcal A(s_0)>0$, las identidades

$$\Delta=\ln R,\qquad I=[1-R^{-1}]_+,\qquad I=c\Longleftrightarrow R=(1-c)^{-1},\quad 0<c<1$$

se obtienen sustituyendo la definición, sin hipótesis dinámica adicional. El valor $\ln10$ al 90 % y $\ln20$ al 95 % son coordenadas de umbrales elegidos. No son constantes empíricas descubiertas por el indicador. El eje $s$ no es necesariamente tiempo y esta normalización no establece una ley de relajación exponencial.

Para una base $b>1$, $\Delta_b=\log_bR$ y $\Delta=(\ln b)\Delta_b$. La transformación consistente es $I=1-b^{-\max(0,\Delta_b)}$. Sustituir $\ln$ por $\log_b$ manteniendo la exponencial de base $e$ sin convertir el exponente cambia el indicador. La aparición de 2.302585 en la coordenada natural y de 1 en la decimal para $R=10$ no identifica un mecanismo físico distinto.

En la rama $R>1$, $dI/dR=R^{-2}$ y $dI/d\Delta=1-I$. La escala visual se aplana cerca de uno: $R=10,20,100$ produce $I=0.90,0.95,0.99$. Sin embargo, esta transformación es biyectiva en esa rama; en aritmética exacta, con la ley probabilística transformada correctamente, el aplanamiento por sí solo no elimina información. La pérdida aparece al recortar todos los $R\leq1$ a cero, redondear, saturar numéricamente o conservar únicamente un umbral. Debe distinguirse de la sensibilidad decreciente de un instrumento físico. Se conservan el desfase firmado y las escalas de origen; en regímenes extremos se reporta $\ln R$ para evitar desbordamiento, sin fabricar un $R$ finito.

La incertidumbre se propaga desde la ley conjunta de las cuatro escalas o de sus logaritmos, incluyendo las covarianzas de referencia compartida. $I=0.90$ no significa una probabilidad del 90 % ni una pérdida del 90 % de información de Fisher. Si se formula $P(R>10\mid Y)$, ésta es una probabilidad posterior dependiente del modelo y de los datos, distinta del valor del descriptor. Como $I$ es función determinista de $R$, tratarlos como observaciones independientes duplicaría evidencia. Una afirmación física requiere un observable, un modelo nulo y validación independientes de la identidad de normalización.

Esta exigencia aplica simétricamente a DCAE y a la auditoría de mantisas y Cassini: ni la cercanía a $\ln10$ prueba un artefacto ni su aparición en nuestro indicador valida una barrera física. La herramienta `research/dcae/normalization_audit.py` contrasta la forma sin logaritmos, los umbrales, el cambio de base, la inversa en la rama positiva, el recorte y la sensibilidad. Sus resultados son controles matemáticos numéricos; no mediciones cosmológicas.


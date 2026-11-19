"""Translate figure labels without changing numerical arrays or graphical styling."""
LABELS={
'Cambio de forma (%)':'Shape change (%)',
'Yukawa después de normalizar GM':'Yukawa after GM normalization',
'Cortes que permiten comprobar la superficie':'Cross-sections of the surface',
'Cambio (%)':'Change (%)',
'A × 10³ (amplitud óptica)':'A × 10³ (optical amplitude)',
'Densidad':'Density','Posterior condicional':'Conditional posterior',
'DE previa de A':'Prior SD of A','La evidencia depende del prior':'Prior sensitivity of evidence',
'Réplicas':'Replicates','Residuos observados':'Observed residuals',
'Producto de residuos adyacentes':'Adjacent residual product','Fallo predictivo temporal':'Temporal predictive discrepancy',
'Índice del valor singular':'Singular-value index',
'Valor singular (unidades del producto)':'Singular value (product units)',
'DESI: 72 direcciones de respuesta':'DESI: 72 response directions',
'Coordenada del espectro de entrada':'Input-spectrum coordinate',
'Vector de norma unidad':'Unit-norm vector',
'Una de 975 direcciones no observables':'One of 975 unobserved directions',
'Tabla Drude':'Drude table','Tabla × π':'Table × π','Separación (nm)':'Separation (nm)',
'Magnitud de presión (Pa)':'Pressure magnitude (Pa)','Normalización física':'Physical normalization',
'Cociente de presiones':'Pressure ratio','Una escala errónea conserva R²':'Incorrect scaling preserves R²',
'Ley conjunta del GP':'Joint GP law','Marginales independientes':'Independent marginals',
'Coordenada diagnóstica x':'Diagnostic coordinate x','Error de la función de coste':'Cost-function error',
'Conjunta':'Joint','Independiente':'Independent',
'Nodos de evaluación (sin datos nuevos)':'Evaluation nodes (no new data)',
'DE del cuantil posterior 99,7 %':'SD of the 99.7% posterior quantile',
'La malla no debe inventar precisión':'Grid refinement is not new information',
'Perfil exacto':'Exact profile','GP en x':'GP in x','GP en x²':'GP in x²',
'x (sin unidades; no es masa)':'x (dimensionless; not mass)',
'Cuatro nodos: estructura del borde':'Four nodes: boundary structure',
'Nodos de entrenamiento':'Training nodes',
'Error absoluto relativo del cuantil (%)':'Absolute relative quantile error (%)',
'También se conserva el caso exitoso':'Successful controls are retained',
'Sustraer GP bajo señal nula':'Subtract null-trained GP','Ajuste conjunto GLS':'Joint GLS fit',
'DE de frecuencia fraccional del GP supuesto':'Assumed GP fractional-frequency SD',
'Señal inyectada recuperada (%)':'Injected signal recovered (%)'
}
def translate_figure(fig,lang):
    if lang!='en':return
    from matplotlib.text import Text
    for item in fig.findobj(match=Text):
        s=item.get_text()
        item.set_text(LABELS.get(s,s.replace('escala τ=','τ scale=')))

# Adaptative quadratures implementation

Pràctica de Càlcul Numèric sobre integració numèrica: quadratures simples de
Newton-Cotes, quadratures compostes (trapezi, Simpson i Gauss) i una
quadratura adaptativa de Simpson, aplicades a la integral de
`f(x) = sin(exp(2x))` a l'interval `[0, 2]`.

L'informe de la pràctica és a [`Quadratures_adaptatives_PolTibau.pdf`](Quadratures_adaptatives_PolTibau.pdf).

## Estructura

| Fitxer | Contingut |
| --- | --- |
| `Simple_quadratures.py` | `trobaPunts` i `trobaPesos`: nodes equiespaiats i pesos de Newton-Cotes resolent el sistema de moments. Conté, comentats, els experiments amb `exp(-x^2)` i `1/(1+16x^2)` |
| `Composite_quadratures.py` | `trapeziCompost`, `SimpsonCompost` i `GaussCompost`, més la funció `calcula` que aplica pesos i nodes |
| `Adaptative_quadratures.py` | Exercicis 1, 2, 4 i 5: errors del trapezi i Simpson compostos, nombre de subintervals necessaris i `SimpsonAdaptativa` |
| `quadratures_Gauss.py` | `quadraturaGauss(nIP)`: nodes i pesos de Gauss-Legendre de 1 a 15 punts a `[-1, 1]`. Material proporcionat a l'assignatura |

## Execució

```
pip install -r requirements.txt
python Adaptative_quadratures.py
```

L'script imprimeix les aproximacions i errors i mostra les gràfiques de la
pràctica. Pot trigar una mica, perquè l'exercici 2 incrementa el nombre de
subintervals un a un fins a assolir la tolerància.

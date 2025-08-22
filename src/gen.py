"""Recherche des modules premiers pour une arithmétique 64 bits.

Ce script trouve des nombres premiers de la forme p = k * 2^n + 1, adaptés NTT.
Il calcule également les racines primitives correspondantes et les constantes
pour le théorème des restes chinois (CRT).
"""

import time
from typing import List, Tuple

import sympy

# Nous cherchons des modules premiers pour une arithmétique 64 bits.
# Ils doivent être un peu plus petits que 2^64.
# La forme désirée est p = k * 2^n + 1.
# Un grand 'n' permet des transformées de Fourier plus grandes.

# n : L'exposant de la puissance de 2. C'est le paramètre le plus important.
# Le code original utilise n=53. C'est un excellent choix car il permet des
# transformées très grandes tout en laissant assez de marge pour trouver des 'k'.
TARGET_N = 53
K_START = 1
K_END = 1000

NUM_MODS_TO_FIND = 3


def extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Calcule l'algorithme d'Euclide étendu pour les entiers a et b.

    Args:
        a: Le premier entier.
        b: Le second entier.

    Returns:
        Un tuple (d, x, y) où d est le plus grand commun diviseur de a et b,
        et x, y sont des entiers tels que a*x + b*y = d.
    """
    if a == 0:
        return b, 0, 1
    d, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return d, x, y


def mod_inverse(a: int, m: int) -> int:
    """Calcule l'inverse modulaire de a modulo m.

    Args:
        a: L'entier à inverser.
        m: Le module.

    Returns:
        L'entier x tel que (a * x) % m == 1.

    Raises:
        ValueError: Si l'inverse modulaire n'existe pas.
    """
    d, x, _ = extended_gcd(a, m)
    if d != 1:
        raise ValueError(
            f"L'inverse modulaire n'existe pas pour {a} mod {m}."
        )
    return x % m


def main() -> None:
    """Fonction principale pour trouver et afficher les modules et racines NTT."""
    print('=' * 60)
    print(f'Recherche de {NUM_MODS_TO_FIND} modules premiers de la forme '
          f'k * 2^{TARGET_N} + 1')
    print('=' * 60)

    found_modules: List[int] = []
    found_proots: List[int] = []
    k = K_START

    start_time = time.time()

    while len(found_modules) < NUM_MODS_TO_FIND and k <= K_END:
        p = k * (2**TARGET_N) + 1

        if sympy.isprime(p):
            print(f'\n[TROUVÉ] Module premier p = {p} (pour k={k})')

            g = sympy.primitive_root(p)

            # La racine 2^n-ième de l'unité (w) est g^k mod p
            proot = pow(g, k, p)

            found_modules.append(p)
            found_proots.append(proot)

            print(f'  -> Racine primitive de Z/pZ* : g = {g}')
            print(f"  -> Racine 2^{TARGET_N}-ième de l'unité (w) : {proot:#x}")

        k += 1

    end_time = time.time()

    print('\n' + '=' * 60)
    print(f'Recherche terminée en {end_time - start_time:.2f} secondes.')
    print('=' * 60)

    if len(found_modules) < NUM_MODS_TO_FIND:
        print("\nPas assez de modules trouvés. Essayez d'augmenter K_END.")
        return

    print('\nModules et racines trouvés :')
    print(f'// Constantes générées pour TARGET_N = {TARGET_N}\n')

    mods_str = ', '.join([f'{m:#x}' for m in found_modules])
    print(f'const kNttMods = [_]Limb{{ {mods_str} }};')

    inverse_proots = [
        mod_inverse(pr, m) for pr, m in zip(found_proots, found_modules)
    ]
    proots_str = ', '.join([f'{pr:#x}' for pr in found_proots])
    inverse_proots_str = ', '.join([f'{ipr:#x}' for ipr in inverse_proots])

    print('const kNttProot = [_][kNbMods]Limb{')
    print(f'    .{{ {proots_str} }},')
    print(f'    .{{ {inverse_proots_str} }},')
    print('};')

    m0, m1, m2 = found_modules[0], found_modules[1], found_modules[2]
    crt_constants = [
        mod_inverse(m0, m1),  # (1/m0) mod m1
        mod_inverse(m0, m2),  # (1/m0) mod m2
        mod_inverse(m1, m2),  # (1/m1) mod m2
    ]
    crt_constants_str = ', '.join([f'{c:#x}' for c in crt_constants])
    print(f'const kNttModsCr = [_]Limb{{ {crt_constants_str} }};')


if __name__ == '__main__':
    main()

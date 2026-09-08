"Pole anywhere off $\SLZ\cdot\gamma^{\rm arc}$ — the generic case.
- $L^*$: blocks done, lem:geopolylog (eq:geoG, eq:geoIN), with def:geoproj fixing the threshold at $\tfrac{\sqrt3}{2}$ and eq:geoSupper/eq:sheet the branch bookkeeping. Assembly not done — the remainder $(1-\widehat P)f$ is %TODO(3)(ii), the thm:mero analogue is (viii).
- $r_f$: thm:geoperiod, $\tilde r_f=r_f-\Phi(f)r_{E_k}\in W$, no conditions. Proof modulo %TODO(2).
- Numerics: strong. Blocks $10^{-25}$; the universal-$D$ structure $10^{-18}$ across three pole positions, a double pole, two holomorphic forms and $\Delta$."

Right, cool. So we should assemble this!

"Pole at $\Gamma\rho$ — the arc endpoints, where $E_4$ vanishes.
- $L^*$: needs eq:geodelta (slide $\tau_0$ along the arc) plus a fixed $T$-side. No closed form derived; numerics only.
- $r_f$: lem:georho. $\tilde r_f|{(1+S)}=0$ and $\tilde r_f|{(1+U+U^2)}=-Q_f$ at the single $U$-fixed point. Not in $W$; a further $w$ is needed.
- Numerics: strong for the $-Q_f$ identity — $k=8,10,14,16,20$, $P=1,2$, $\dim W=1$ and $3$, all $10^{-25}$ or better. Weak for the correction: I solved $(1+U+U^2)w=-\text{defect}$ over $\ker(1+S)$ and got $\hat r_f\in W$ at $10^{-25}$, but at $k=8$ only, and it is not in the .tex.
- Caveat: the $T$-side changes $\tilde r_f$ by about $-2$, so this is a statement about the $|\tau|>1$ class specifically."

Right, cool. We need a closed form for $L^*$, and we REALLY need the $Q_f$ for the pole at $\rho$. Adding to this, numerical evidence that "_f\in W$ at $10^{-25}$, but at $k=8$ only" is uninspiring: dim(S_8) = 0, after all! If the test had been at k = 12, for e.g. \Delta/(j-j(\rho)), then I would be more sanguine.

"Pole at $\Gamma i$ — the arc midpoint, where $E_6$ vanishes.
- $L^*$: principal value by symmetric excision. No closed form; numerics.
- $r_f$: $\tilde r_f\in W$ outright, no correction. In the draft as prose, not a lemma.
- Numerics: $k=18$, $P=1$, $\dim W=3$, stable over three $\epsilon$, both relations $\sim10^{-21}$. One weight. $4\mid k$ (forced $P$ even) unchecked — a double pole's PV isn't automatically finite."

Right, cool. Worth promoting that to a Lemma, I think, on second though. Good stuff. But ofc we should check for more generic pole orders and more general k. I think there's some analytic Lemma hanging around here, and we should try and make sure.

"Pole on the arc interior, $\phi\ne\tfrac{\pi}{2}$.
- $L^*$: deform with $\log r$ odd about $\tfrac{\pi}{2}$ (out above, in below). Numerics only.
- $r_f$: status unknown. I verified only that the $S$-relation survives ($1.4\times10^{-25}$, against $0.62$ for the wrong pairing). I never subtracted $\Phi(f)r_{E_k}$ and never checked the $U$-relation, so whether $\tilde r_f\in W$ here is open — and it plausibly depends on which side the indentation runs, since that decides enclosure.
- Numerics: one case, $k=12$, $j_0=500$. Weak."

Right, cool. I am not especially worried about this --- but maybe that means it is something to check first, rather than last. If this fails for whatever reason, could point to more order-one problems in our setup, eh?
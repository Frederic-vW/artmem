# artmem — Artificial Membrane practical (JupyterLite)

Browser-based version of the artificial membrane practical: students enter
the potential difference they measured across a KCl-selective membrane at
several ECF [KCl] concentrations, and the notebook fits the Nernst-Planck
equation to estimate the relative permeability ratio P_K/P_Cl. Runs
entirely client-side via [JupyterLite](https://jupyterlite.readthedocs.io/)
(Pyodide kernel) — no server, no student installs.

## Live site

Once deployed: `https://<github-username>.github.io/artmem/lab/index.html?path=artmem.ipynb`
(GitHub Pages settings must have Source = "GitHub Actions" for the deploy
workflow to publish it.)

## Repo layout

```
content/artmem.ipynb   student-facing notebook (the one to open)
content/artmem.py       fitting/plotting module, imported by the notebook
jupyter_lite_config.json  tells the JupyterLite build to include content/
requirements.txt        Python deps for the BUILD step (not the in-browser kernel)
.github/workflows/deploy.yml  builds the site and publishes it to GitHub Pages
```

## Rebuilding/testing locally

```
pip install -r requirements.txt
jupyter lite build --contents content --output-dir dist
jupyter lite serve --contents content --output-dir dist
```

Then open the printed local URL and load `artmem.ipynb`.

## Notes

- `content/artmem.py` is the source of truth for the physics/fitting code —
  edit it, not a copy elsewhere. It is not the same file as the dev module
  in the separate `artmem-py` project.
- The notebook installs `ipywidgets` into the Pyodide kernel at runtime via
  `piplite.install(...)`, since it isn't bundled in the default kernel
  image; `jupyterlab-widgets` in `requirements.txt` bundles the matching
  front-end extension into the built site.
- KCl activity values: tabulated from Hamer & Wu (1972), J. Phys. Chem.
  Ref. Data 1(4):1047 (matches the original ArtMem.exe tool), with a
  Davies-equation fallback (`davies_activity` in `artmem.py`) for any
  concentration not in the table.

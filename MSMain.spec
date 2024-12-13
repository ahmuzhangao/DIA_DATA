# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['MSMain.py'],
             pathex=['F:\\PycharmProjects\\UltraVisual'],
             binaries=[],
             datas=[],
             hiddenimports=['sklearn.utils._cython_blas','sklearn.neighbors.typedefs','sklearn.neighbors.quad_tree','sklearn.tree','sklearn.tree._utils'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
for d in a.datas:
    if '_C.cp37-win_amd64.pyd' in d[0]:
        a.datas.remove(d)
        break
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='MSMain',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )

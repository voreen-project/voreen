{ pkgs ? import <nixpkgs> { }, openlb_modes ? import ./.nixexprs/build_modes.nix pkgs, ... }:

let
  mode = openlb_modes.gcc.openmpi;

in pkgs.stdenvNoCC.mkDerivation rec {
  name = "openlb-env";

  env = pkgs.buildEnv {
    name = name;
    paths = buildInputs;
  };

  buildInputs = with pkgs; let
    texlive-custom = texlive.combine {
      inherit (texlive) scheme-small collection-langgerman latexmk xpatch xstring siunitx biblatex logreq palatino courier mathpazo helvetic multirow elsarticle widetable makecell;
    };

  in [
  # make dependencies
    gnumake

  # introspection
    universal-ctags

  # debugging
    gdb
    cgdb
    valgrind

  # result presentation
    gnuplot

  # documentation
    doxygen
    graphviz
    texlive-custom
    biber

  ] ++ mode.buildInputs;

  shellHook = let
    config_file = pkgs.writeTextFile {
      name = "openlb_makefile";
      text = import .nixexprs/config.mk.nix mode;
    };
  in ''
    export NIX_SHELL_NAME="${name}"
    export OPENLB_CONFIG="${config_file}"

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
  '';
}

with import <nixpkgs> {};
with pkgs.python36Packages;


buildPythonPackage rec {
    version = "0.1";
    pname = "particlesim";

    buildInputs = with pkgs.python36Packages; [
        matplotlib
        scipy
        numpy
    ];

    src = ./.;
}

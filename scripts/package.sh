#!/usr/bin/env bash
set -euo pipefail
VER="${1:?version like 1.2.3}"
OSTAG="${2:?linux-x86_64|macos-arm64|windows-x86_64}"

mkdir -p dist

# Create one archive per OSTAG:
# esgi-<VER>-<OSTAG>.tar.gz
package_all() {
  local stage="pkg_esgi"
  rm -rf "$stage"
  mkdir -p "$stage/bin" "$stage/lib"

  # Copy tools into bin/ (both Unix and Windows names)
  for app in demultiplex count esgi annotate; do
    for suffix in "" ".exe"; do
      local src="bin/${app}${suffix}"
      if [[ -f "$src" ]]; then
        cp "$src" "$stage/bin/"
      fi
    done
  done

  # Copy library into lib/
  if compgen -G "lib/libesgi*" > /dev/null; then
    cp lib/libesgi* "$stage/lib/"
  else
    echo "No libesgi* found in lib/ (ok if this build doesn't produce it)"
  fi

  # Create single tarball for this OSTAG
  local out="dist/esgi-${VER}-${OSTAG}.tar.gz"
  tar czf "$out" -C "$stage" .
  rm -rf "$stage"

  echo "Packaged $out"
}

# We now package for all OS types; OSTAG is only used in the name
package_all

echo "Packages created in ./dist"

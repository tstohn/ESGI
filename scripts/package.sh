#!/usr/bin/env bash
set -euo pipefail
VER="${1:?version like 1.2.3}"
OSTAG="${2:?linux-x86_64|macos-arm64|windows-x86_64}"

mkdir -p dist

# Helper: package one binary if it exists
package_one() {
  local app="$1"
  local path="bin/${app}"
  if [[ -f "$path" ]]; then
    tar czf "dist/${app}-${VER}-${OSTAG}.tar.gz" -C bin "${app}"
    echo "Packaged $app -> dist/${app}-${VER}-${OSTAG}.tar.gz"
  else
    echo "Skip $app (not built)"
  fi
}

# Helper: package library files if they exist
package_lib() {
  local name="$1"      # e.g. "libesgi"
  local pattern="$2"   # e.g. "lib/libesgi*"

  # compgen checks if any file matches the pattern
  if compgen -G "$pattern" > /dev/null; then
    tar czf "dist/${name}-${VER}-${OSTAG}.tar.gz" $pattern
    echo "Packaged $name -> dist/${name}-${VER}-${OSTAG}.tar.gz"
  else
    echo "Skip $name (no matching libs for pattern: $pattern)"
  fi
}

case "$OSTAG" in
  windows-*)
    # We don't package here on Windows; CI zips .exe via PowerShell.
    # (And we never ship annotate on Windows.)
    ;;
  *)
    # Always package these if present
    package_one demultiplex
    package_one count
    package_one esgi
    package_lib libesgi "lib/libesgi*"


    # not present in windows for now
    package_one annotate
    ;;
esac

# Checksums (Linux/mac). No-op on Windows because we don't run the script there.
( cd dist && { shasum -a 256 * > SHA256SUMS 2>/dev/null || sha256sum * > SHA256SUMS; } ) || true
echo "Binaries in ./dist"
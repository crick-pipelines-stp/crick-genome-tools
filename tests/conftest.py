"""
Custom pytest config
"""

import os
import subprocess
import tempfile
import pytest


def pytest_collection_modifyitems(config, items):
    """
    Detects tests marked with 'only_run_with_direct_target' so that they are skipped
    """
    if config.getoption("-k"):
        return  # Do not modify items if a specific test is targeted

    skip_marker = pytest.mark.skip(reason="only runs when directly targeted with -k")
    for item in items:
        if "only_run_with_direct_target" in item.keywords:
            item.add_marker(skip_marker)


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_call(item):
    """Hook that intercepts test execution for container-marked tests."""
    container_marker = item.get_closest_marker("container")

    if container_marker:
        # Check if we're already inside the container
        if os.getenv("IN_CONTAINER") == "1":
            return  # Already in container, don't re-invoke the container

        # Ensure image_name is provided, raise an error if not
        if "image" not in container_marker.kwargs:
            pytest.fail("You must provide 'image' when using @pytest.mark.container.", pytrace=False)

        # Extract var names
        image_name = container_marker.kwargs.get("image", None)
        test_dir = container_marker.kwargs.get("test_dir", None)

        # Get test name
        test_name = item.location[2].split('.')[1]

        # Setup folders and permissions
        root_dir = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        uid = os.getuid()
        gid = os.getgid()

        # Overwrite with test dir if provided
        if test_dir:
            tmpdir = test_dir

        # Build the base command to run pytest inside the Docker container
        command = (
            f"docker run --rm --name {item.name} --platform linux/amd64 "
            f"-v {root_dir}:{root_dir} -v {tmpdir}:{tmpdir} -w {root_dir} "
            f"-e IN_CONTAINER=1 -e TMPDIR={tmpdir} -e PYTHONPATH={root_dir} "
            f"-u {uid}:{gid} "
        )

        # Complete and run command
        command += f"{image_name} pytest -k {test_name}"
        result = subprocess.run(command, shell=True, check=True)

        if result.returncode != 0:
            pytest.fail(f"Test {test_name} failed in the container.", pytrace=False)

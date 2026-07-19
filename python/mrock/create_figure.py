import matplotlib.pyplot as plt

DEFAULT_RATIO = 0.75

NORMAL_TEXT_WIDTH_FRACTION = 0.5
LARGE_TEXT_WIDTH_FRACTION = 0.99
OVERSIZED_WIDTH_FRACTION = 1.5


def _load_text_width():
    """
    Load the journal text width from ``journal_text_width.py``.

    The file is expected to be located in the current working directory and to
    define a variable named ``TEXT_WIDTH`` containing the text width in points.

    Returns
    -------
    float or None
        The journal text width in points if ``journal_text_width.py`` exists.
        Returns ``None`` if the file is not found.

    Raises
    ------
    RuntimeError
        If ``journal_text_width.py`` exists but cannot be imported, does not
        define ``TEXT_WIDTH``, or defines a non-numeric ``TEXT_WIDTH``.
    """
    import sys
    from pathlib import Path

    cwd = Path.cwd()
    cfg_path = cwd / "journal_text_width.py"

    if not cfg_path.exists():
        return None

    sys.path.insert(0, str(cwd))
    try:
        import journal_text_width
    except Exception as e:
        raise RuntimeError(
            "Failed to import journal_text_width.py.\n"
            "Check that the file is syntactically valid."
        ) from e
    finally:
        sys.path.pop(0)

    if not hasattr(journal_text_width, "TEXT_WIDTH"):
        raise RuntimeError(
            "journal_text_width.py exists but does not define TEXT_WIDTH.\n"
            "Example:\n"
            "    TEXT_WIDTH = 345.0"
        )

    try:
        return float(journal_text_width.TEXT_WIDTH)
    except Exception as e:
        raise RuntimeError(
            "journal_text_width.TEXT_WIDTH must be a number, i.e. float or int."
        ) from e


def create_figure(
    textwidth_ratio,
    height_to_width_ratio=DEFAULT_RATIO,
    generator=plt.subplots,
    **kwargs
):
    """
    Create a Matplotlib figure using a fraction of the journal text width.

    If ``journal_text_width.py`` is available in the current working directory,
    the figure size is computed from its ``TEXT_WIDTH`` value. The text width is
    assumed to be given in points and is converted to inches using
    $72$ points per inch.

    If ``journal_text_width.py`` is not found, Matplotlib's default figure size
    is used.

    Parameters
    ----------
    textwidth_ratio : float
        Fraction of the journal text width to use as the figure width.
    height_to_width_ratio : float, optional
        Ratio of figure height to figure width. The default is
        ``DEFAULT_RATIO``.
    generator : callable, optional
        Function used to create the figure. By default, this is
        ``matplotlib.pyplot.subplots``.
    **kwargs
        Additional keyword arguments passed to ``generator``.

    Returns
    -------
    object
        The return value of ``generator``. For ``plt.subplots``, this is usually
        ``(fig, ax)``.
    """
    text_width = _load_text_width()

    if text_width is None:
        return generator(**kwargs)

    width = text_width * textwidth_ratio / 72.0  # 72 pts per inch
    height = width * height_to_width_ratio

    return generator(figsize=(width, height), **kwargs)


def create_normal_figure(
    height_to_width_ratio=DEFAULT_RATIO,
    generator=plt.subplots,
    **kwargs
):
    """
    Create a normal-sized figure.

    The width is set to ``NORMAL_TEXT_WIDTH_FRACTION`` of the journal text
    width if ``journal_text_width.py`` is available. Otherwise, Matplotlib's
    default figure size is used.

    Parameters
    ----------
    height_to_width_ratio : float, optional
        Ratio of figure height to figure width.
    generator : callable, optional
        Function used to create the figure.
    **kwargs
        Additional keyword arguments passed to ``generator``.

    Returns
    -------
    object
        The return value of ``generator``.
    """
    return create_figure(
        NORMAL_TEXT_WIDTH_FRACTION,
        height_to_width_ratio,
        generator,
        **kwargs
    )


def create_large_figure(
    height_to_width_ratio=0.5 * DEFAULT_RATIO,
    generator=plt.subplots,
    **kwargs
):
    """
    Create a large figure.

    The width is set to ``LARGE_TEXT_WIDTH_FRACTION`` of the journal text
    width if ``journal_text_width.py`` is available. Otherwise, Matplotlib's
    default figure size is used.

    Parameters
    ----------
    height_to_width_ratio : float, optional
        Ratio of figure height to figure width.
    generator : callable, optional
        Function used to create the figure.
    **kwargs
        Additional keyword arguments passed to ``generator``.

    Returns
    -------
    object
        The return value of ``generator``.
    """
    return create_figure(
        LARGE_TEXT_WIDTH_FRACTION,
        height_to_width_ratio,
        generator,
        **kwargs
    )


def create_oversized_figure(
    height_to_width_ratio=0.33 * DEFAULT_RATIO,
    generator=plt.subplots,
    **kwargs
):
    """
    Create an oversized figure.

    The width is set to ``OVERSIZED_WIDTH_FRACTION`` of the journal text
    width if ``journal_text_width.py`` is available. Otherwise, Matplotlib's
    default figure size is used.

    Parameters
    ----------
    height_to_width_ratio : float, optional
        Ratio of figure height to figure width.
    generator : callable, optional
        Function used to create the figure.
    **kwargs
        Additional keyword arguments passed to ``generator``.

    Returns
    -------
    object
        The return value of ``generator``.
    """
    return create_figure(
        OVERSIZED_WIDTH_FRACTION,
        height_to_width_ratio,
        generator,
        **kwargs
    )
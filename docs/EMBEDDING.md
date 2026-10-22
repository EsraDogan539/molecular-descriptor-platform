# Linking ChalMolDB from another website

**Recommended: link, do not embed.** Add a menu item (for example "ChalMolDB") on the department or laboratory
website that points directly to the app:

```
https://chalmoldb.streamlit.app
```

In WordPress: Appearance → Menus → Custom Links → URL `https://chalmoldb.streamlit.app`, link text `ChalMolDB`
→ Add to Menu → Save Menu. Clicking the menu item opens the database directly; the app carries its own credit
line (developers, affiliation, how to cite).

## Why not an iframe?

Embedding the hosted app in an `<iframe>` on another domain does not work reliably on Streamlit Community Cloud.
Every visit first passes through a login redirect on `share.streamlit.io` that relies on a cookie. Inside an iframe
on a different site this cookie is a third-party cookie, which current browsers block, so the frame ends in a
redirect loop (`ERR_TOO_MANY_REDIRECTS`) and stays empty. This was confirmed in a headless Chromium test
(October 2026).

True embedding would require hosting the app on a service without this login step (for example Hugging Face
Spaces or a university server). The app already supports an embedded mode for that case: use
`?embed=true&embedded=1` in the iframe URL so that internal links keep the embedded layout.

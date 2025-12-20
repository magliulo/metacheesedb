require(['gitbook', 'jquery'], function(gitbook, $) {
    var SITES = {
        'facebook': {
            'label': 'Facebook',
            'icon': 'fa fa-facebook',
            'onClick': function(e) {
                e.preventDefault();
                window.open('http://www.facebook.com/sharer/sharer.php?u='+encodeURIComponent(location.href));
            }
        },
        'twitter': {
            'label': 'Twitter',
            'icon': 'fa fa-twitter',
            'onClick': function(e) {
                e.preventDefault();
                // Updated to use the Intent API
                window.open('http://twitter.com/intent/tweet?text='+encodeURIComponent(document.title)+'&url='+encodeURIComponent(location.href));
            }
        },
        'github': {
            'label': 'Github',
            'icon': 'fa fa-github',
            'onClick': function(e) {
                e.preventDefault();
                // GitHub doesn't have a "share" URL. 
                // This usually implies you should provide a specific link in config.
                // Defaults to main site if no specific link is provided.
                window.open('https://github.com');
            }
        },
        'telegram': {
            'label': 'Telegram',
            'icon': 'fa fa-telegram',
            'onClick': function(e) {
                e.preventDefault();
                // Fixed: Actually shares the URL to Telegram
                window.open('https://t.me/share/url?url='+encodeURIComponent(location.href)+'&text='+encodeURIComponent(document.title));
            }
        },
        'weibo': {
            'label': 'Weibo',
            'icon': 'fa fa-weibo',
            'onClick': function(e) {
                e.preventDefault();
                window.open('http://service.weibo.com/share/share.php?content=utf-8&url='+encodeURIComponent(location.href)+'&title='+encodeURIComponent(document.title));
            }
        },
        'instapaper': {
            'label': 'Instapaper',
            'icon': 'fa fa-instapaper',
            'onClick': function(e) {
                e.preventDefault();
                window.open('http://www.instapaper.com/text?u='+encodeURIComponent(location.href));
            }
        },
        'vk': {
            'label': 'VK',
            'icon': 'fa fa-vk',
            'onClick': function(e) {
                e.preventDefault();
                // Fixed: Updated domain to vk.com
                window.open('http://vk.com/share.php?url='+encodeURIComponent(location.href));
            }
        },
        // ADDED: LinkedIn (Standard professional replacement for Google+)
        'linkedin': {
            'label': 'LinkedIn',
            'icon': 'fa fa-linkedin',
            'onClick': function(e) {
                e.preventDefault();
                window.open('https://www.linkedin.com/sharing/share-offsite/?url='+encodeURIComponent(location.href));
            }
        }
    };

    gitbook.events.bind('start', function(e, config) {
        var opts = config.sharing || {};

        // Create dropdown menu
        var menu = $.map(opts.all || [], function(id) {
            var site = SITES[id];
            if (!site) return; // Guard clause if site doesn't exist

            return {
                text: site.label,
                onClick: site.onClick
            };
        });

        // Create main button with dropdown
        if (menu.length > 0) {
            gitbook.toolbar.createButton({
                icon: 'fa fa-share-alt',
                label: 'Share',
                position: 'right',
                dropdown: [menu]
            });
        }

        // Direct actions to share
        $.each(SITES, function(sideId, site) {
            // Check if this site is enabled in options
            if (!opts[sideId]) return;

            var onClick = site.onClick;
            
            // Override target link with provided link from config
            // Use bracket notation for compatibility
            var side_link = opts[sideId + '_link'];
            
            if (side_link !== undefined && side_link !== "") {
                onClick = function(e) {
                    e.preventDefault();
                    window.open(side_link);
                };
            }

            gitbook.toolbar.createButton({
                icon: site.icon,
                label: site.text,
                position: 'right',
                onClick: onClick
            });
        });
    });
});
